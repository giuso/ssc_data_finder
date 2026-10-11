"""GUI page listing everything downloadable for one GEO series (opened via "?downloads=GSExxxx")."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from ssc_data_finder.geo import geo_url
from ssc_data_finder.geo_details import (
    HOW_TO_LOAD,
    SeriesDetails,
    fetch_details,
    list_folder,
    series_folder,
    size_in_bytes,
)


def downloads_url(accession: str) -> str:
    """Link into this app that opens the downloads page for `accession` (new tab from a LinkColumn)."""
    return f"?downloads={accession}#{accession}"


@st.cache_data(ttl=24 * 3600, show_spinner=False)
def _details(accession: str) -> SeriesDetails:
    return fetch_details(accession)


@st.cache_data(ttl=24 * 3600, show_spinner=False)
def _folder(url: str) -> dict[str, str]:
    return list_folder(url)


def render(accession: str) -> None:
    accession = accession.strip().upper()
    if not accession.startswith("GSE") or not accession[3:].isdigit():
        st.error(f"'{accession}' is not a GEO series accession (GSE followed by digits).")
        return

    with st.spinner(f"Reading {accession} from GEO…"):
        try:
            d = _details(accession)
            folder = series_folder(accession)
            suppl_sizes = _folder(folder + "suppl/")
            metadata = {kind: _folder(folder + f"{kind}/") for kind in ("matrix", "soft", "miniml")}
        except Exception as exc:
            st.error(f"Could not read {accession} from GEO: {exc}")
            return

    st.title(f"Downloads · {accession}")
    st.caption(d.title)
    st.markdown(
        " · ".join(
            [
                f"[GEO page]({geo_url(accession)})",
                f"[GEO download page](https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession})",
                f"[FTP folder]({folder})",
                f"[SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc={accession})",
                "[← New search](?)",
            ]
        )
    )

    sample_files = [(s, f) for s in d.samples for f in s.supplementary_files]
    total = sum(size_in_bytes(suppl_sizes.get(f.name)) for f in d.supplementary_files)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Processed files", len(d.supplementary_files))
    c2.metric("Total size", _human(total) if total else "—")
    c3.metric("Samples", len(d.samples))
    c4.metric("Per-sample files", len(sample_files))

    # --- Series-level processed data -------------------------------------------
    st.subheader("Processed data")
    if d.supplementary_files:
        st.dataframe(
            pd.DataFrame(
                {
                    "File": [f.url for f in d.supplementary_files],
                    "Format": [f.format for f in d.supplementary_files],
                    "Size": [suppl_sizes.get(f.name, "") for f in d.supplementary_files],
                    "How to load": [HOW_TO_LOAD.get(f.format, "") for f in d.supplementary_files],
                }
            ),
            hide_index=True,
            width="stretch",
            column_config={
                "File": st.column_config.LinkColumn(display_text=r"/([^/]+)$", help="Click to download", width="large"),
                "How to load": st.column_config.TextColumn(width="large"),
            },
        )
    else:
        st.write("No series-level files; see the per-sample files and raw reads below.")

    # --- Per-sample files ------------------------------------------------------
    st.subheader("Per-sample files")
    if sample_files:
        if any(f.format == "Raw archive (tar)" for f in d.supplementary_files):
            st.caption(f"The same files are bundled in {accession}_RAW.tar above.")
        formats = sorted({f.format for _, f in sample_files})
        chosen = st.multiselect("Format", formats, placeholder="All formats")
        rows = [(s, f) for s, f in sample_files if not chosen or f.format in chosen]
        st.dataframe(
            pd.DataFrame(
                {
                    "Sample": [geo_url(s.accession) for s, _ in rows],
                    "Sample title": [s.title for s, _ in rows],
                    "Condition": [s.value_of("disease") for s, _ in rows],
                    "File": [f.url for _, f in rows],
                    "Format": [f.format for _, f in rows],
                }
            ),
            hide_index=True,
            width="stretch",
            column_config={
                "Sample": st.column_config.LinkColumn(display_text=r"acc=(GSM\d+)"),
                "File": st.column_config.LinkColumn(display_text=r"/([^/]+)$", help="Click to download", width="large"),
            },
        )
    else:
        st.write("No per-sample files.")

    # --- Metadata --------------------------------------------------------------
    st.subheader("Sample metadata")
    labels = {
        "matrix": "Series Matrix: sample characteristics as a table",
        "soft": "SOFT family file: full series + sample records (text)",
        "miniml": "MINiML: the same records as XML",
    }
    meta_rows = [
        (folder + f"{kind}/{name}", labels[kind], size) for kind, files in metadata.items() for name, size in files.items()
    ]
    if meta_rows:
        st.dataframe(
            pd.DataFrame(meta_rows, columns=["File", "Contents", "Size"]),
            hide_index=True,
            width="stretch",
            column_config={"File": st.column_config.LinkColumn(display_text=r"/([^/]+)$", width="large")},
        )

    # --- Raw reads -------------------------------------------------------------
    st.subheader("Raw sequencing reads")
    st.markdown(
        f"FASTQ/BAM files are stored in SRA, not GEO. Open the "
        f"[SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc={accession}) to see every run "
        f"with its size, then download the accession list and fetch the runs with the SRA Toolkit."
    )
    st.code("prefetch --option-file SRR_Acc_List.txt\nfasterq-dump --split-files SRR_ACCESSION", language="bash")

    # --- Command line ----------------------------------------------------------
    st.subheader("Download from the command line")
    lines = [
        f"# Everything in the processed-data folder (saved under {accession}/suppl/)",
        f'wget -r -np -nH --cut-dirs=3 -R "index.html*" {folder}suppl/',
    ]
    if d.supplementary_files:
        lines += ["", "# Or file by file"] + [f"wget {f.url}" for f in d.supplementary_files]
    st.code("\n".join(lines), language="bash")


def _human(size: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1000 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{size:.0f} B"
        size /= 1000
    return ""
