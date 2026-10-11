"""Common organism names from MeSH tags, GEO taxa and abstract text, so papers and datasets filter alike."""
from __future__ import annotations

import re

NOT_STATED = "Not stated"

# Scientific name (GEO taxon) -> common name.
TAXA = {
    "Homo sapiens": "Human",
    "Mus musculus": "Mouse",
    "Rattus norvegicus": "Rat",
    "Danio rerio": "Zebrafish",
    "Sus scrofa": "Pig",
    "Macaca mulatta": "Macaque",
    "Macaca fascicularis": "Macaque",
    "Canis lupus familiaris": "Dog",
    "Gallus gallus": "Chicken",
    "Drosophila melanogaster": "Fruit fly",
    "Bos taurus": "Cow",
    "Oryctolagus cuniculus": "Rabbit",
}

# MeSH descriptor (prefix, so "Mice, Inbred C57BL" counts) -> common name.
MESH = {
    "Humans": "Human",
    "Mice": "Mouse",
    "Rats": "Rat",
    "Zebrafish": "Zebrafish",
    "Swine": "Pig",
    "Macaca": "Macaque",
    "Dogs": "Dog",
    "Chickens": "Chicken",
    "Drosophila": "Fruit fly",
    "Cattle": "Cow",
    "Rabbits": "Rabbit",
}

# Fallback for papers without MeSH tags (recent, not yet indexed): words in the abstract.
TEXT = {
    "Human": r"\b(humans?|patients?|healthy (donors?|controls?|volunteers?|individuals?)|participants?)\b",
    "Mouse": r"\b(mouse|mice|murine)\b",
    "Rat": r"\brats?\b",
    "Zebrafish": r"\bzebrafish\b",
    "Pig": r"\b(pigs?|porcine|swine)\b",
    "Macaque": r"\b(macaques?|non-human primates?|nonhuman primates?)\b",
}


def ordered(names) -> list[str]:
    """Distinct names, common lab organisms first (Human, Mouse, Rat, ...), then alphabetical."""
    order = list(dict.fromkeys(TAXA.values()))
    return sorted(dict.fromkeys(names), key=lambda n: (n not in order, order.index(n) if n in order else 0, n))


def from_taxon(taxon: str | None) -> list[str]:
    """"Homo sapiens; Mus musculus" -> ["Human", "Mouse"]."""
    names = [t.strip() for t in re.split(r";|,", taxon or "") if t.strip()]
    return list(dict.fromkeys(TAXA.get(n, n) for n in names))


def from_mesh(descriptors: list[str]) -> list[str]:
    found = [common for d in descriptors for prefix, common in MESH.items() if d.startswith(prefix)]
    return list(dict.fromkeys(found))


def from_text(text: str) -> list[str]:
    return [name for name, pattern in TEXT.items() if re.search(pattern, text or "", re.IGNORECASE)]
