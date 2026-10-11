"""Fetch and parse the full GEO record of a series and its samples (GEO "SOFT" text format)."""
from __future__ import annotations

import pickle
import re
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import requests

ACC_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
CACHE_DIR = Path.home() / ".cache" / "ssc_data_finder" / "geo"
CACHE_MAX_AGE = 7 * 24 * 3600  # GEO records rarely change once public
CACHE_VERSION = "v2"  # bump when SeriesDetails/Sample gain fields, so old pickles are refetched
FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/geo"

# Sample characteristic keys are free text; these map common spellings onto one field.
FIELD_ALIASES = {
    "tissue": ["tissue", "tissue type", "organ", "tissue source", "anatomical site", "biopsy site", "site"],
    "disease": [
        "disease state", "disease", "condition", "diagnosis", "disease status", "health status", "health_status",
        "subject status", "sample condition", "status", "group", "phenotype",
    ],
    "cell_type": ["cell type", "cell population", "cell subset", "sorted population", "tissue/cell type"],
    "sex": ["sex", "gender"],
    "age": ["age", "age (years)", "age (yrs)", "age at biopsy", "age at sampling"],
    "treatment": ["treatment", "agent", "drug", "stimulation"],
}

SAMPLE_ID_RE = re.compile(r"[A-Za-z_-]*\d+[\w-]*")  # "Hama04", "SC266": a sample name, not a tissue
HEALTHY_RE = re.compile(r"^(healthy|normal|control|hc|nc)\b|\b(healthy|normal) (control|donor|subject|volunteer)s?\b", re.I)

# Checked against sample protocols/processing text, first match wins.
TECHNOLOGIES = [
    ("10x Visium (spatial)", r"visium|spatial gene expression"),
    ("10x Multiome", r"multiome"),
    ("10x Chromium", r"10\s?x|chromium|cell ?ranger"),
    ("Smart-seq", r"smart-?seq"),
    ("BD Rhapsody", r"rhapsody"),
    ("Drop-seq", r"drop-?seq"),
    ("inDrop", r"indrop"),
    ("Seq-Well", r"seq-?well"),
    ("Parse Biosciences", r"parse bio|split-?seq|evercode"),
    ("CEL-seq", r"cel-?seq"),
    ("MARS-seq", r"mars-?seq"),
]

# Supplementary file name -> what you get when you download it.
FILE_FORMATS = [
    ("Seurat object", r"\.(rds|h5seurat|robj|rdata)(\.gz)?$"),
    ("AnnData (h5ad)", r"\.h5ad(\.gz)?$"),
    ("10x HDF5", r"\.h5(\.gz)?$"),
    ("10x matrix (mtx)", r"(\.mtx|barcodes\.tsv|features\.tsv|genes\.tsv)(\.gz)?$"),
    ("Loom", r"\.loom(\.gz)?$"),
    ("Counts table", r"\.(csv|tsv|txt|xlsx?)(\.gz)?$"),
    ("Raw archive (tar)", r"_raw\.tar$|\.tar(\.gz)?$"),
]


@dataclass
class Sample:
    accession: str
    title: str
    source: str | None
    organism: str | None
    characteristics: dict[str, str]
    instrument: str | None
    library_source: str | None
    library_strategy: str | None
    technology: str | None
    supplementary_files: list[SupplementaryFile] = field(default_factory=list)

    def value_of(self, name: str) -> str | None:
        """Harmonised value for a FIELD_ALIASES name, e.g. sample.value_of("tissue")."""
        value = next((self.characteristics[a] for a in FIELD_ALIASES[name] if a in self.characteristics), None)
        if value is None and name == "tissue" and self.source and not SAMPLE_ID_RE.fullmatch(self.source):
            value = self.source
        return tidy_value(value, name) if value else None


@dataclass
class SupplementaryFile:
    name: str
    url: str
    format: str


@dataclass
class SeriesDetails:
    accession: str
    title: str
    summary: str
    overall_design: str
    series_type: list[str]
    submission_date: str | None
    public_date: str | None
    pubmed_ids: list[str]
    contributors: list[str]
    contact_name: str | None
    institute: str | None
    city: str | None
    country: str | None
    platforms: list[str]
    organisms: list[str]
    supplementary_files: list[SupplementaryFile]
    relations: dict[str, str]
    samples: list[Sample] = field(default_factory=list)

    @property
    def submission_year(self) -> str | None:
        match = re.search(r"\d{4}", self.submission_date or "")
        return match.group(0) if match else None

    @property
    def technologies(self) -> list[str]:
        return _distinct(s.technology for s in self.samples)

    @property
    def instruments(self) -> list[str]:
        return _distinct(s.instrument for s in self.samples)

    @property
    def file_formats(self) -> list[str]:
        return _distinct(f.format for f in self.supplementary_files)

    def values(self, name: str) -> Counter:
        """Sample counts per harmonised value, e.g. details.values("disease")."""
        return Counter(v for s in self.samples if (v := s.value_of(name)))


def tidy_value(value: str, name: str | None = None) -> str:
    """Make spelling variants count as one: "skin" -> "Skin"; healthy/normal controls -> "Healthy control"."""
    value = re.sub(r"\s+", " ", value).strip(" _")
    if name == "disease" and HEALTHY_RE.search(value):
        return "Healthy control"
    return value[:1].upper() + value[1:]


def parse_soft(text: str) -> list[dict[str, list[str]]]:
    """Split SOFT text into records ("^SERIES", "^SAMPLE", ...), each {"!Key": [values]}."""
    records: list[dict[str, list[str]]] = []
    for line in text.splitlines():
        if line.startswith("^"):
            records.append({})
        elif line.startswith("!") and records and " = " in line:
            key, value = line.split(" = ", 1)
            records[-1].setdefault(key, []).append(value.strip())
    return records


def cached_details(accession: str) -> SeriesDetails | None:
    """The on-disk copy of a series fetched within the last week, if any."""
    cache_file = CACHE_DIR / f"{accession}.{CACHE_VERSION}.pickle"
    if cache_file.exists() and time.time() - cache_file.stat().st_mtime < CACHE_MAX_AGE:
        try:
            return pickle.loads(cache_file.read_bytes())
        except Exception:  # stale format after an upgrade: refetch
            return None
    return None


def fetch_details(accession: str, use_cache: bool = True) -> SeriesDetails:
    """Full GEO record of a series and all its samples, cached on disk for a week."""
    if use_cache and (details := cached_details(accession)):
        return details
    cache_file = CACHE_DIR / f"{accession}.{CACHE_VERSION}.pickle"

    series = parse_soft(_fetch_soft(accession, "self"))[0]
    details = _series_from_record(accession, series)
    details.samples = [_sample_from_record(r) for r in parse_soft(_fetch_soft(accession, "gsm"))]

    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file.write_bytes(pickle.dumps(details))
    return details


def fetch_many(accessions: list[str], workers: int = 3, on_done=None) -> dict[str, SeriesDetails | Exception]:
    """fetch_details for several series in parallel; failures are returned, not raised.

    `on_done(accession, result)` is called as each one finishes (e.g. to drive a progress bar).
    """
    out: dict[str, SeriesDetails | Exception] = {}

    def task(acc: str):
        try:
            return acc, fetch_details(acc)
        except Exception as exc:  # one broken record shouldn't sink the rest
            return acc, exc

    # NCBI asks for at most ~3 requests/second without an API key.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for acc, result in pool.map(task, accessions):
            out[acc] = result
            if on_done:
                on_done(acc, result)
    return out


# How to open each downloadable format.
HOW_TO_LOAD = {
    "Seurat object": "R: readRDS() for .rds; SeuratDisk::LoadH5Seurat() for .h5Seurat",
    "AnnData (h5ad)": "Python: scanpy.read_h5ad()",
    "10x HDF5": "Python: scanpy.read_10x_h5() · R: Seurat::Read10X_h5()",
    "10x matrix (mtx)": "Python: scanpy.read_10x_mtx() · R: Seurat::Read10X() (needs matrix + barcodes + features)",
    "Loom": "Python: scanpy.read_loom() / loompy",
    "Counts table": "Python: pandas.read_csv() · R: read.csv() / data.table::fread()",
    "Raw archive (tar)": "All per-sample files bundled; extract with: tar -xf <file>",
    "Other": "",
}


def series_folder(accession: str) -> str:
    """GEO's FTP folder for a series: GSE264508 -> .../geo/series/GSE264nnn/GSE264508/."""
    number = accession[3:]
    return f"{FTP_ROOT}/series/GSE{number[:-3]}nnn/{accession}/"


def list_folder(url: str) -> dict[str, str]:
    """File name -> human-readable size ("1.9G") from a GEO FTP folder listing (served over https)."""
    response = requests.get(url, timeout=30)
    if response.status_code == 404:
        return {}
    response.raise_for_status()
    pattern = r'<a href="([^"?/][^"]*)">[^<]*</a>\s+\d{4}-\d\d-\d\d \d\d:\d\d\s+([\d.]+[KMGT]?)'
    return dict(re.findall(pattern, response.text))


def size_in_bytes(size: str | None) -> float:
    """"1.9G" -> 1.9e9 (for totals and sorting)."""
    if not size:
        return 0.0
    units = {"K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
    return float(size.rstrip("KMGT")) * units.get(size[-1], 1)


def https_url(url: str) -> str:
    """GEO lists supplementary files as ftp:// links; the same paths are served over https."""
    return re.sub(r"^ftp://", "https://", url)


def _fetch_soft(accession: str, target: str, retries: int = 4) -> str:
    params = {"acc": accession, "targ": target, "form": "text", "view": "brief"}
    for attempt in range(retries):
        response = requests.get(ACC_URL, params=params, timeout=60)
        if response.status_code in (429, 500, 502, 503, 504) and attempt < retries - 1:
            time.sleep(1 + 2 * attempt)
            continue
        response.raise_for_status()
        response.encoding = "utf-8"  # GEO sends UTF-8 without declaring it; requests would assume Latin-1
        return response.text
    raise AssertionError("unreachable")


def _series_from_record(accession: str, r: dict[str, list[str]]) -> SeriesDetails:
    def one(key: str) -> str | None:
        return r.get(key, [None])[0]

    relations = {}
    for rel in r.get("!Series_relation", []):
        if ": " in rel:
            kind, target = rel.split(": ", 1)
            relations.setdefault(kind, target)

    files = []
    for url in r.get("!Series_supplementary_file", []):
        if url.upper() == "NONE":
            continue
        name = url.rsplit("/", 1)[-1]
        files.append(SupplementaryFile(name=name, url=https_url(url), format=_classify_file(name)))

    return SeriesDetails(
        accession=accession,
        title=one("!Series_title") or "",
        summary=" ".join(r.get("!Series_summary", [])),
        overall_design=" ".join(r.get("!Series_overall_design", [])),
        series_type=r.get("!Series_type", []),
        submission_date=one("!Series_submission_date"),
        public_date=(one("!Series_status") or "").replace("Public on ", "") or None,
        pubmed_ids=r.get("!Series_pubmed_id", []),
        contributors=[_person(c) for c in r.get("!Series_contributor", [])],
        contact_name=_person(one("!Series_contact_name") or "") or None,
        institute=one("!Series_contact_institute"),
        city=one("!Series_contact_city"),
        country=one("!Series_contact_country"),
        platforms=r.get("!Series_platform_id", []),
        organisms=_distinct(r.get("!Series_sample_organism", []) or r.get("!Series_platform_organism", [])),
        supplementary_files=files,
        relations=relations,
    )


def _sample_from_record(r: dict[str, list[str]]) -> Sample:
    def one(key: str) -> str | None:
        return r.get(key, [None])[0]

    characteristics = {}
    for item in r.get("!Sample_characteristics_ch1", []):
        key, sep, value = item.partition(": ")
        if sep and value:
            characteristics.setdefault(key.strip().lower(), value.strip())

    protocol_text = " ".join(
        r.get("!Sample_extract_protocol_ch1", [])
        + r.get("!Sample_data_processing", [])
        + r.get("!Sample_description", [])
        + list(characteristics.values())
    )
    return Sample(
        accession=one("!Sample_geo_accession") or "",
        title=one("!Sample_title") or "",
        source=one("!Sample_source_name_ch1"),
        organism=one("!Sample_organism_ch1"),
        characteristics=characteristics,
        instrument=one("!Sample_instrument_model"),
        library_source=one("!Sample_library_source"),
        library_strategy=one("!Sample_library_strategy"),
        supplementary_files=[
            SupplementaryFile(name=url.rsplit("/", 1)[-1], url=https_url(url), format=_classify_file(url))
            for key, urls in r.items()
            if key.startswith("!Sample_supplementary_file")
            for url in urls
            if url.upper() != "NONE"
        ],
        technology=next(
            (name for name, pattern in TECHNOLOGIES if re.search(pattern, protocol_text, re.IGNORECASE)), None
        ),
    )


def _classify_file(name: str) -> str:
    return next((label for label, pattern in FILE_FORMATS if re.search(pattern, name, re.IGNORECASE)), "Other")


def _person(raw: str) -> str:
    """GEO writes names as "First,Middle,Last" ("Kathryn,,Torok")."""
    return " ".join(p for p in raw.split(",") if p)


def _distinct(values) -> list[str]:
    return list(dict.fromkeys(v for v in values if v))
