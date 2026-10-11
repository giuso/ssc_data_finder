"""Search Europe PMC for papers using single-cell RNA-seq techniques."""
from __future__ import annotations

import html
import re
from dataclasses import dataclass, field

import requests

from ssc_data_finder import organisms
from ssc_data_finder.query import quote_topic, split_authors

EUROPEPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

SCRNA_SEQ_SYNONYMS = [
    "single-cell RNA sequencing",
    "single cell RNA sequencing",
    "single-cell RNA-seq",
    "scRNA-seq",
    "single-cell transcriptomics",
]

GEO_ACCESSION_RE = re.compile(r"\bGSE\d+\b")

# Affiliation segments are tried in this order; the first keyword that matches
# a comma-separated segment names the institution ("Dept of X, Y University, City").
INSTITUTION_KEYWORDS = [
    "university",
    "universit",  # Università, Universität, Université, Universidad
    "institute",
    "institut",
    "hospital",
    "hôpital",
    "klinik",
    "clinic",
    "centre",
    "center",
    "college",
    "school of medicine",
    "foundation",
    "laboratory",
]
EMAIL_RE = re.compile(r"\S+@\S+")
COUNTRY_ALIASES = {
    "people's republic of china": "China",
    "pr china": "China",
    "p.r. china": "China",
    "p. r. china": "China",
    "united states": "USA",
    "united states of america": "USA",
    "u.s.a": "USA",
    "us": "USA",
    "united kingdom": "UK",
    "republic of korea": "Korea",
    "south korea": "Korea",
    "the netherlands": "Netherlands",
}
DEPARTMENT_RE = re.compile(r"(department|dept|division|section|program|unit)\b", re.IGNORECASE)
SHORT_INSTITUTION_RE = re.compile(
    r"[A-Z][\w'&-]*\s(?:University|School of Medicine)|University of\s[A-Z][\w'-]*"
)


@dataclass
class Author:
    name: str
    institutions: list[str] = field(default_factory=list)


@dataclass
class Paper:
    pmid: str | None
    doi: str | None
    title: str
    journal: str | None
    pub_year: str | None
    abstract: str
    geo_accessions: list[str] = field(default_factory=list)
    authors: list[Author] = field(default_factory=list)
    mesh_organisms: list[str] = field(default_factory=list)  # from MeSH tags; empty if not yet indexed


def build_query(topic: str | None = None, author: str | None = None) -> str:
    synonyms = " OR ".join(f'"{s}"' for s in SCRNA_SEQ_SYNONYMS)
    clauses = []
    if topic and topic.strip():
        quoted = quote_topic(topic)
        if quoted != topic.strip():
            # Plain topic: match it in title/abstract/keywords only - Europe PMC otherwise also
            # searches full text, pulling in papers that just cite the topic.
            quoted = f"TITLE:{quoted} OR ABSTRACT:{quoted} OR KW:{quoted}"
        clauses.append(f"({quoted})")
    authors = split_authors(author)
    if authors:
        clauses.append("(" + " OR ".join(f'AUTH:"{a}"' for a in authors) + ")")
    if not clauses:
        raise ValueError("Give a topic, an author, or both.")
    return " AND ".join(clauses) + f" AND ({synonyms})"


def normalize_institution(affiliation: str) -> str | None:
    """Reduce one free-text affiliation to its institution, e.g. "Central South University, China"."""
    text = EMAIL_RE.sub("", affiliation)
    text = re.sub(r"electronic address:?", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\([^)]*\)", "", text)
    text = re.sub(r"\s+", " ", text)
    segments = [seg.strip(" .") for seg in re.split(r",|\.\s", text) if seg.strip(" .")]
    if not segments:
        return None
    country = segments[-1] if len(segments) > 1 and not re.search(r"[\d.@/]", segments[-1]) else None
    if country:
        country = COUNTRY_ALIASES.get(country.lower(), country)
    for keyword in INSTITUTION_KEYWORDS:
        for seg in segments:
            if keyword not in seg.lower():
                continue
            seg = re.sub(r"^(the|and)\s+", "", seg, flags=re.IGNORECASE)
            if len(segments) == 1 or DEPARTMENT_RE.match(seg):
                # No commas, e.g. "Department of Medicine Yale University New Haven Connecticut USA":
                # keep just the short "Yale University" / "University of X" phrase.
                match = SHORT_INSTITUTION_RE.search(seg)
                if match:
                    return match.group(0)
            return f"{seg}, {country}" if country and country != seg else seg
    return None


def parse_authors(record: dict) -> list[Author]:
    authors = []
    for a in record.get("authorList", {}).get("author", []):
        if a.get("firstName") and a.get("lastName"):
            name = f"{a['firstName']} {a['lastName']}"
        else:
            name = a.get("fullName") or a.get("collectiveName")
        if not name:
            continue
        affiliations = a.get("authorAffiliationDetailsList", {}).get("authorAffiliation", [])
        institutions = []
        for aff in affiliations:
            # One entry can hold several affiliations separated by ";".
            for part in aff.get("affiliation", "").split(";"):
                inst = normalize_institution(part)
                if inst and inst not in institutions:
                    institutions.append(inst)
        authors.append(Author(name=name, institutions=institutions))
    return authors


def clean_text(text: str) -> str:
    """Unescape HTML entities and drop inline tags such as <sup>/<i> from titles."""
    return re.sub(r"<[^>]+>", "", html.unescape(text or "")).strip()


def search_papers(topic: str | None = None, max_results: int = 25, author: str | None = None) -> list[Paper]:
    """Search Europe PMC for scRNA-seq papers matching `topic` and/or written by `author`."""
    params = {
        "query": build_query(topic, author),
        "format": "json",
        "resultType": "core",
        "pageSize": min(max_results, 1000),
    }
    response = requests.get(EUROPEPMC_SEARCH_URL, params=params, timeout=30)
    response.raise_for_status()
    results = response.json().get("resultList", {}).get("result", [])

    papers = []
    for r in results[:max_results]:
        abstract = r.get("abstractText", "") or ""
        papers.append(
            Paper(
                pmid=r.get("pmid"),
                doi=r.get("doi"),
                title=clean_text(r.get("title", "")),
                journal=r.get("journalInfo", {}).get("journal", {}).get("title") or r.get("journalTitle"),
                pub_year=r.get("pubYear"),
                abstract=abstract,
                geo_accessions=sorted(set(GEO_ACCESSION_RE.findall(abstract))),
                authors=parse_authors(r),
                mesh_organisms=organisms.from_mesh(
                    [h.get("descriptorName", "") for h in r.get("meshHeadingList", {}).get("meshHeading", [])]
                ),
            )
        )
    return papers
