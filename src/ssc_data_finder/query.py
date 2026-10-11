"""Helpers shared by the Europe PMC and GEO query builders."""
from __future__ import annotations

import re


def quote_topic(topic: str) -> str:
    """Search the topic as an exact phrase, unless it is already a hand-written query
    (contains quotes or AND/OR/NOT), so "systemic sclerosis" doesn't match "systemic" + "sclerosis"."""
    topic = topic.strip()
    if '"' in topic or re.search(r"\b(AND|OR|NOT)\b", topic):
        return topic
    return f'"{topic}"'


def split_authors(author: str | None) -> list[str]:
    """"John Varga; Varga J" -> ["John Varga", "Varga J"]."""
    return [a.strip() for a in re.split(r"[;\n]", author or "") if a.strip()]


def last_initials(name: str) -> str:
    """Normalise to the PubMed/GEO "Last FM" form: "John Varga" -> "Varga J"; "Varga J" stays."""
    parts = name.replace(",", " ").split()
    if len(parts) < 2 or re.fullmatch(r"[A-Z]{1,3}", parts[-1]):
        return " ".join(parts)
    if "," in name:  # "Varga, John"
        last, firsts = parts[0], parts[1:]
    else:
        last, firsts = parts[-1], parts[:-1]
    # Hyphenated first names give one initial per part: "Pei-Suen Tsou" -> "Tsou PS".
    return f"{last} {''.join(p[0].upper() for f in firsts for p in f.split('-') if p)}"


def author_matches(name: str, query: str) -> bool:
    """True if a full author name ("John Varga") matches a searched name.

    "Varga J" matches any J. Varga; "John Varga" requires the first name to match too.
    """
    full, wanted = last_initials(name).lower().split(), last_initials(query).lower().split()
    if not full or not wanted or full[0] != wanted[0]:
        return False
    if len(wanted) == 1 or len(full) == 1:
        return True
    if full[1][0] != wanted[1][0]:
        return False
    query_parts = query.replace(",", " ").split()
    if last_initials(query) == " ".join(query_parts):  # query was already "Varga J"
        return True
    first = query_parts[1] if "," in query else query_parts[0]
    return name.lower().startswith(first.lower())
