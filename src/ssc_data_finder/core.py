"""Core search + cross-reference logic, shared by the CLI, the GUI and the notebook."""
from __future__ import annotations

import re
from difflib import SequenceMatcher

from ssc_data_finder import europepmc, geo, organisms

# GEO series are often titled after their paper but never get its PMID filled in.
TITLE_MATCH_THRESHOLD = 0.9


def _normalize_title(title: str) -> str:
    return re.sub(r"[^a-z0-9 ]", "", re.sub(r"[-–—_/]", " ", title.lower())).strip()


def titles_match(a: str, b: str) -> bool:
    a, b = _normalize_title(a), _normalize_title(b)
    return bool(a and b) and SequenceMatcher(None, a, b).ratio() >= TITLE_MATCH_THRESHOLD


def cross_reference(
    papers, geo_series, pmid_links: dict[str, list[str]] | None = None, all_series: list | None = None
) -> list[dict]:
    """One row per paper. A GEO series counts as linked to a paper when:

    - NCBI links it to the paper's PMID (`pmid_links`: PMID -> accessions), or
    - the GEO record lists the paper's PMID, or
    - the abstract mentions an accession that the GEO search also found, or
    - the GEO series has (nearly) the same title as the paper.

    `organisms` comes from the paper's MeSH tags plus its linked datasets' organisms
    (`all_series` supplies those), falling back to words in the abstract.
    """
    organism_of = {s.accession: s.organisms for s in (all_series or geo_series)}
    pmid_links = pmid_links or {}
    searched = {s.accession for s in geo_series}
    by_pmid: dict[str, set[str]] = {}
    for s in geo_series:
        for pmid in s.pubmed_ids:
            by_pmid.setdefault(pmid, set()).add(s.accession)

    rows = []
    for p in papers:
        linked = set(p.geo_accessions) & searched
        if p.pmid:
            linked |= set(pmid_links.get(p.pmid, [])) | by_pmid.get(p.pmid, set())
        linked |= {s.accession for s in geo_series if titles_match(s.title, p.title)}
        species = organisms.ordered(p.mesh_organisms + [o for a in linked for o in organism_of.get(a, [])])
        species = species or organisms.from_text(f"{p.title} {p.abstract}")
        rows.append(
            {
                "title": p.title,
                "journal": p.journal,
                "pub_year": p.pub_year,
                "pmid": p.pmid,
                "doi": p.doi,
                "authors": [{"name": a.name, "institutions": a.institutions} for a in p.authors],
                "institutions": sorted({i for a in p.authors for i in a.institutions}),
                "mentioned_geo_accessions": p.geo_accessions,
                "linked_geo_accessions": sorted(linked),
                "has_known_dataset": bool(linked),
                "organisms": species,
            }
        )
    return rows


def dataset_rows(series: list[geo.GeoSeries], linked_accessions: set[str]) -> list[dict]:
    """One row per GEO series, flagging whether its paper is among the Europe PMC results."""
    return [
        {
            "accession": s.accession,
            "url": s.url,
            "title": s.title,
            "organism": s.organism,
            "organisms": s.organisms,
            "n_samples": s.n_samples,
            "pubmed_ids": s.pubmed_ids,
            "paper_in_results": s.accession in linked_accessions,
        }
        for s in sorted(series, key=lambda s: s.accession, reverse=True)
    ]


def search_all(topic: str | None = None, max_results: int = 25, author: str | None = None) -> dict:
    """Search Europe PMC + GEO and cross-reference both ways.

    Returns {"papers": [...], "datasets": [...]}: the papers (see `cross_reference`) and every GEO
    series found, either by the GEO search itself or because NCBI links it to one of the papers.
    """
    papers = europepmc.search_papers(topic, max_results=max_results, author=author)
    pmids = [p.pmid for p in papers if p.pmid]

    search_uids = geo.search_uids(topic, max_results=max_results, author=author)
    uid_links = geo.linked_uids(pmids)
    series = geo.fetch_series(sorted(set(search_uids) | {u for uids in uid_links.values() for u in uids}))

    pmid_links = {
        pmid: [series[u].accession for u in uids if u in series] for pmid, uids in uid_links.items()
    }
    searched = [series[u] for u in search_uids if u in series]
    rows = cross_reference(papers, searched, pmid_links, all_series=list(series.values()))
    linked = {a for row in rows for a in row["linked_geo_accessions"]}
    return {"papers": rows, "datasets": dataset_rows(list(series.values()), linked)}


def run_search(topic: str | None = None, max_results: int = 25, author: str | None = None) -> list[dict]:
    """Search Europe PMC + GEO by `topic` and/or `author`; return the cross-referenced paper rows.

    Several authors can be given separated by ";" (matched with OR).
    """
    return search_all(topic, max_results=max_results, author=author)["papers"]
