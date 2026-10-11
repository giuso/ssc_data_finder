"""Search NCBI GEO (gds database) for deposited single-cell RNA-seq datasets."""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import requests

from ssc_data_finder import organisms
from ssc_data_finder.query import last_initials, quote_topic, split_authors

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
GEO_ACCESSION_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}"

SCRNA_SEQ_SYNONYMS = [
    "single-cell RNA sequencing",
    "single cell RNA-seq",
    "scRNA-seq",
    "single-cell transcriptomics",
]


@dataclass
class GeoSeries:
    accession: str
    title: str
    summary: str
    organism: str | None
    dataset_type: str | None
    n_samples: int | None = None
    pubmed_ids: list[str] = field(default_factory=list)

    @property
    def url(self) -> str:
        return geo_url(self.accession)

    @property
    def organisms(self) -> list[str]:
        """Common names, e.g. ["Human", "Mouse"]."""
        return organisms.from_taxon(self.organism)


def geo_url(accession: str) -> str:
    """GEO accession page, e.g. GSE264508 -> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264508."""
    return GEO_ACCESSION_URL.format(accession)


def build_term(topic: str | None = None, author: str | None = None) -> str:
    synonyms = " OR ".join(f'"{s}"' for s in SCRNA_SEQ_SYNONYMS)
    clauses = []
    if topic and topic.strip():
        clauses.append(f"({quote_topic(topic)})")
    authors = split_authors(author)
    if authors:
        # GEO only matches the "Last F" form ("Varga J"), not "John Varga".
        clauses.append("(" + " OR ".join(f'"{last_initials(a)}"[AUTH]' for a in authors) + ")")
    if not clauses:
        raise ValueError("Give a topic, an author, or both.")
    return " AND ".join(clauses) + f" AND ({synonyms})"


def search_uids(
    topic: str | None = None, max_results: int = 25, api_key: str | None = None, author: str | None = None
) -> list[str]:
    """GEO (gds) UIDs of single-cell RNA-seq entries matching `topic` and/or submitted by `author`."""
    params = {"db": "gds", "term": build_term(topic, author), "retmode": "json", "retmax": max_results}
    result = _eutils("get", ESEARCH_URL, params, api_key)
    return result.get("esearchresult", {}).get("idlist", [])


def linked_uids(pmids: list[str], api_key: str | None = None) -> dict[str, list[str]]:
    """GEO (gds) UIDs that NCBI links to each PubMed ID - the datasets deposited with that paper."""
    if not pmids:
        return {}
    # POST with one `id` per PMID so NCBI returns a separate linkset per paper.
    params = {"dbfrom": "pubmed", "db": "gds", "retmode": "json", "id": pmids}
    result = _eutils("post", ELINK_URL, params, api_key)
    links = {}
    for linkset in result.get("linksets", []):
        uids = [u for db in linkset.get("linksetdbs", []) for u in db.get("links", [])]
        if uids:
            links[linkset["ids"][0]] = uids
    return links


def fetch_series(uids: list[str], api_key: str | None = None) -> dict[str, GeoSeries]:
    """Summaries for GEO UIDs, keyed by UID; only series (GSE) are kept."""
    if not uids:
        return {}
    params = {"db": "gds", "id": ",".join(uids), "retmode": "json"}
    result = _eutils("post", ESUMMARY_URL, params, api_key).get("result", {})

    series = {}
    for uid in result.get("uids", []):
        doc = result[uid]
        accession = doc.get("accession") or f"GSE{doc.get('gse', '')}"
        if not accession.startswith("GSE"):
            continue
        series[uid] = GeoSeries(
            accession=accession,
            title=doc.get("title", ""),
            summary=doc.get("summary", ""),
            organism=doc.get("taxon"),
            dataset_type=doc.get("gdstype"),
            n_samples=doc.get("n_samples"),
            pubmed_ids=[str(p) for p in doc.get("pubmedids", [])],
        )
    return series


def search_datasets(
    topic: str | None = None, max_results: int = 25, api_key: str | None = None, author: str | None = None
) -> list[GeoSeries]:
    """Search GEO for single-cell RNA-seq series matching `topic` and/or submitted by `author`."""
    uids = search_uids(topic, max_results=max_results, api_key=api_key, author=author)
    return list(fetch_series(uids, api_key=api_key).values())


def _eutils(method: str, url: str, params: dict, api_key: str | None, retries: int = 3) -> dict:
    """Call an E-utilities endpoint, backing off on NCBI's 3-requests/second limit (HTTP 429)."""
    if api_key:
        params = {**params, "api_key": api_key}
    for attempt in range(retries):
        if method == "post":
            response = requests.post(url, data=params, timeout=30)
        else:
            response = requests.get(url, params=params, timeout=30)
        if response.status_code != 429 or attempt == retries - 1:
            response.raise_for_status()
            return response.json()
        time.sleep(1 + attempt)
    raise AssertionError("unreachable")
