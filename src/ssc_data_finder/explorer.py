"""GUI: explore the full GEO records (series + samples) of the datasets found by a search."""
from __future__ import annotations

from collections import Counter
from dataclasses import replace

import altair as alt
import pandas as pd
import streamlit as st

from ssc_data_finder.downloads import downloads_url
from ssc_data_finder.geo import geo_url
from ssc_data_finder import organisms
from ssc_data_finder.geo_details import Sample, SeriesDetails, cached_details, fetch_many

# Single-series bars: categorical slot 1 (blue), its own step in dark mode.
BAR = {"light": "#2a78d6", "dark": "#3987e5"}
TEXT_SECONDARY = {"light": "#52514e", "dark": "#c3c2b7"}
TOP_N = 12
DOWNLOAD_FORMATS = {"Seurat object", "AnnData (h5ad)", "10x HDF5", "10x matrix (mtx)", "Loom", "Counts table"}


def render(datasets: list[dict], mode: str, organism_filter: set[str] | None = None) -> None:
    """`organism_filter`: common names to keep (None = all), from the page-wide filter row."""
    accessions = [d["accession"] for d in datasets]
    if not accessions:
        st.info("No GEO datasets for this search.")
        return

    cache: dict[str, SeriesDetails] = st.session_state.setdefault("geo_details", {})
    for acc in accessions:
        if acc not in cache and (details := cached_details(acc)):
            cache[acc] = details
    missing = [a for a in accessions if a not in cache]
    if missing:
        st.caption(
            f"Reads the full GEO page of each dataset: study design, submitting lab, every sample's "
            f"characteristics, sequencing technology and downloadable files. "
            f"{len(missing)} dataset(s) to fetch (about 2 s each the first time; cached on disk for a week)."
        )
        if not st.button(f"Fetch details for {len(missing)} GEO dataset(s)", type="primary"):
            if len(missing) == len(accessions):
                return
            st.caption(f"Showing the {len(accessions) - len(missing)} dataset(s) already fetched.")
        else:
            _fetch(missing, cache)

    details = [cache[a] for a in accessions if a in cache]
    if details:
        details = _filter_samples(details, organism_filter)
        if not details:
            st.info("No samples match these filters.")
            return
        _overview(details, mode)
        st.divider()
        _single_dataset(details, mode)


def _fetch(accessions: list[str], cache: dict[str, SeriesDetails]) -> None:
    progress = st.progress(0.0, text="Fetching GEO records…")
    done, failed = 0, {}

    def on_done(acc, result):
        nonlocal done
        done += 1
        progress.progress(done / len(accessions), text=f"Fetched {acc} ({done}/{len(accessions)})")

    for acc, result in fetch_many(accessions, on_done=on_done).items():
        if isinstance(result, Exception):
            failed[acc] = result
        else:
            cache[acc] = result
    progress.empty()
    if failed:
        st.warning(
            "Could not fetch " + ", ".join(failed) + " (GEO did not answer); press the button again to retry."
        )


def _filter_samples(details: list[SeriesDetails], organism_filter: set[str] | None) -> list[SeriesDetails]:
    """Sample-level filter row (tissue, condition, technology); keeps datasets with samples left."""

    def organism_ok(s: Sample) -> bool:
        return organism_filter is None or bool(set(organisms.from_taxon(s.organism) or ["Not stated"]) & organism_filter)

    pool = [s for d in details for s in d.samples if organism_ok(s)]
    st.markdown("**Filter samples**")
    cols = st.columns(3)
    chosen = {}
    for col, (name, label) in zip(cols, [("tissue", "Tissue"), ("disease", "Disease / condition"), ("technology", None)]):
        values = Counter(_sample_value(s, name) for s in pool)
        with col:
            chosen[name] = set(
                st.multiselect(
                    label or "Technology",
                    sorted(values, key=lambda v: -values[v]),
                    format_func=lambda v, values=values: f"{v} ({values[v]})",
                    placeholder="All",
                    key=f"explorer_{name}",
                )
            )

    def keep(s: Sample) -> bool:
        return organism_ok(s) and all(not want or _sample_value(s, name) in want for name, want in chosen.items())

    filtered = [replace(d, samples=[s for s in d.samples if keep(s)]) for d in details]
    return [d for d in filtered if d.samples or (not d.samples and not any(chosen.values()) and organism_filter is None)]


def _sample_value(s: Sample, name: str) -> str:
    value = s.technology if name == "technology" else s.value_of(name)
    return value or "Not stated"


# --- Across all datasets -------------------------------------------------------


def _overview(details: list[SeriesDetails], mode: str) -> None:
    samples = sum(len(d.samples) for d in details)
    processed = sum(bool(set(d.file_formats) & DOWNLOAD_FORMATS) for d in details)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Datasets", len(details))
    c2.metric("Samples", samples)
    c3.metric("With processed data to download", processed, help="Seurat, AnnData, 10x matrix/HDF5 or counts table")
    c4.metric("Submitting labs", len({d.institute for d in details if d.institute}))

    def count(name):
        return sum((d.values(name) for d in details), Counter())

    def per_dataset(attr):
        return Counter(v for d in details for v in (getattr(d, attr) or ["Not stated"]))

    charts = [
        ("Samples by tissue", count("tissue"), "samples"),
        ("Samples by disease / condition", count("disease"), "samples"),
        ("Datasets by technology", per_dataset("technologies"), "datasets"),
        ("Datasets by sequencing instrument", per_dataset("instruments"), "datasets"),
        ("Datasets by downloadable file format", per_dataset("file_formats"), "datasets"),
        ("Datasets by submitting lab", Counter(_lab(d) for d in details), "datasets"),
    ]
    for row in range(0, len(charts), 2):
        cols = st.columns(2)
        for col, (title, counts, unit) in zip(cols, charts[row : row + 2]):
            with col:
                st.markdown(f"**{title}**")
                _hbar(counts, unit, mode)

    left, right = st.columns(2)
    with left:
        st.markdown("**Datasets by submission year**")
        _year_bars(Counter(d.submission_year for d in details if d.submission_year), mode)
    with right:
        st.markdown("**Samples by sex**")
        _hbar(count("sex"), "samples", mode)

    st.markdown("**Datasets**")
    st.dataframe(
        pd.DataFrame(
            {
                "Accession": [geo_url(d.accession) for d in details],
                "Title": [d.title for d in details],
                "Samples": [len(d.samples) for d in details],
                "Organism": [", ".join(d.organisms) for d in details],
                "Technology": [", ".join(d.technologies) for d in details],
                "Tissues": [", ".join(d.values("tissue")) for d in details],
                "Conditions": [", ".join(d.values("disease")) for d in details],
                "Files": [", ".join(d.file_formats) for d in details],
                "Lab": [_lab(d) for d in details],
                "Submitted": [d.submission_date for d in details],
            }
        ),
        hide_index=True,
        width="stretch",
        column_config={
            "Accession": st.column_config.LinkColumn(display_text=r"acc=(GSE\d+)"),
            "Title": st.column_config.TextColumn(width="large"),
        },
    )

    samples_df = pd.DataFrame(
        [
            {
                "Dataset": d.accession,
                "Sample": s.accession,
                "Title": s.title,
                "Tissue": s.value_of("tissue"),
                "Disease / condition": s.value_of("disease"),
                "Cell type": s.value_of("cell_type"),
                "Sex": s.value_of("sex"),
                "Age": s.value_of("age"),
                "Treatment": s.value_of("treatment"),
                "Organism": s.organism,
                "Technology": s.technology,
                "Instrument": s.instrument,
            }
            for d in details
            for s in d.samples
        ]
    )
    with st.expander(f"All {len(samples_df)} samples (harmonised characteristics)"):
        st.dataframe(samples_df, hide_index=True, width="stretch")
        st.download_button(
            "Download samples CSV", samples_df.to_csv(index=False).encode(), "geo_samples.csv", "text/csv"
        )


# --- One dataset ---------------------------------------------------------------


def _single_dataset(details: list[SeriesDetails], mode: str) -> None:
    st.markdown("#### Explore one dataset")
    by_acc = {d.accession: d for d in details}
    acc = st.selectbox(
        "Dataset",
        list(by_acc),
        format_func=lambda a: f"{a} · {by_acc[a].title[:110]}",
        label_visibility="collapsed",
    )
    d = by_acc[acc]

    st.markdown(f"##### [{d.accession}]({geo_url(d.accession)}) · {d.title}")
    facts = {
        "Samples": len(d.samples),
        "Organism": ", ".join(d.organisms) or "—",
        "Technology": ", ".join(d.technologies) or "Not stated",
        "Submitted": d.submission_date or "—",
        "Public since": d.public_date or "—",
    }
    # Plain text rather than st.metric: values like "10x Visium (spatial)" are too long for metric tiles.
    for col, (label, value) in zip(st.columns(len(facts)), facts.items()):
        col.caption(label)
        col.markdown(f"**{value}**")

    links = [f"[GEO]({geo_url(d.accession)})", f"[Downloads]({downloads_url(d.accession)})"]
    links += [f"[PubMed {p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)" for p in d.pubmed_ids]
    links += [f"[{kind}]({url})" for kind, url in d.relations.items() if url.startswith("http")]
    links.append(f"[SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc={d.accession})")
    st.markdown(" · ".join(links))

    left, right = st.columns([3, 2])
    with left:
        st.markdown("**Summary**")
        st.write(d.summary or "—")
        if d.overall_design:
            st.markdown("**Overall design**")
            st.write(d.overall_design)
    with right:
        st.markdown("**Submitted by**")
        place = ", ".join(x for x in (d.institute, d.city, d.country) if x)
        st.write(f"{d.contact_name or '—'}  \n{place}")
        if d.contributors:
            st.markdown("**Contributors**")
            st.write(", ".join(d.contributors))
        st.markdown("**Instruments**")
        st.write(", ".join(d.instruments) or "—")

    # One small chart per sample characteristic that actually splits the samples into groups.
    keys = Counter(k for s in d.samples for k in s.characteristics)
    splits = []
    for key in keys:
        values = Counter(s.characteristics.get(key) for s in d.samples if s.characteristics.get(key))
        if 2 <= len(values) <= 15 and len(values) < len(d.samples):
            splits.append((key, values))
    if splits:
        st.markdown("**How the samples split**")
        for row in range(0, len(splits), 3):
            for col, (key, values) in zip(st.columns(3), splits[row : row + 3]):
                with col:
                    st.caption(key.capitalize())
                    _hbar(values, "samples", mode)

    st.markdown("**Samples**")
    sample_keys = list(keys)
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "Sample": geo_url(s.accession),
                    "Title": s.title,
                    "Source": s.source,
                    **{k.capitalize(): s.characteristics.get(k) for k in sample_keys},
                    "Technology": s.technology,
                    "Instrument": s.instrument,
                }
                for s in d.samples
            ]
        ),
        hide_index=True,
        width="stretch",
        column_config={"Sample": st.column_config.LinkColumn(display_text=r"acc=(GSM\d+)")},
    )

    st.markdown("**Supplementary files**")
    if d.supplementary_files:
        st.dataframe(
            pd.DataFrame(
                {"File": [f.url for f in d.supplementary_files], "Format": [f.format for f in d.supplementary_files]}
            ),
            hide_index=True,
            width="stretch",
            column_config={"File": st.column_config.LinkColumn(display_text=r"/([^/]+)$", help="Click to download")},
        )
    else:
        st.write("None listed.")


# --- Charts --------------------------------------------------------------------


def _hbar(counts: Counter, unit: str, mode: str) -> None:
    if not counts:
        st.caption("Not stated in GEO.")
        return
    items = counts.most_common(TOP_N)
    rest = counts.most_common()[TOP_N:]
    if rest:
        items.append((f"Other ({len(rest)} values)", sum(n for _, n in rest)))
    df = pd.DataFrame(items, columns=["label", "count"])
    order = list(df["label"])

    base = alt.Chart(df).encode(
        y=alt.Y("label:N", sort=order, title=None, axis=alt.Axis(labelLimit=220, ticks=False, domain=False)),
        x=alt.X("count:Q", title=None, axis=None),
        tooltip=[alt.Tooltip("label:N", title="Value"), alt.Tooltip("count:Q", title=unit.capitalize())],
    )
    bars = base.mark_bar(color=BAR[mode], cornerRadiusEnd=4)
    labels = base.mark_text(align="left", dx=4, color=TEXT_SECONDARY[mode]).encode(text="count:Q")
    chart = (bars + labels).properties(height=max(60, 24 * len(df))).configure_scale(bandPaddingInner=0.3)
    st.altair_chart(chart.configure_view(stroke=None), use_container_width=True)


def _year_bars(counts: Counter, mode: str) -> None:
    if not counts:
        st.caption("Not stated in GEO.")
        return
    years = range(int(min(counts)), int(max(counts)) + 1)  # keep empty years so gaps read as gaps
    df = pd.DataFrame({"year": [str(y) for y in years], "count": [counts.get(str(y), 0) for y in years]})
    base = alt.Chart(df).encode(
        x=alt.X("year:O", title=None, axis=alt.Axis(labelAngle=0, ticks=False)),
        y=alt.Y("count:Q", title="Datasets", axis=alt.Axis(tickMinStep=1, grid=True, gridOpacity=0.3)),
        tooltip=[alt.Tooltip("year:O", title="Year"), alt.Tooltip("count:Q", title="Datasets")],
    )
    chart = base.mark_bar(color=BAR[mode], cornerRadiusEnd=4).properties(height=220)
    st.altair_chart(chart.configure_scale(bandPaddingInner=0.3).configure_view(stroke=None), use_container_width=True)


def _lab(d: SeriesDetails) -> str:
    return ", ".join(x for x in (d.institute, d.country) if x) or "Not stated"
