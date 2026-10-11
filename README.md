# ssc_data_finder

Find published papers that use single-cell RNA-seq on a given topic (or by a
given author), cross-reference them against datasets deposited in NCBI GEO, and
explore the research centers, authors and datasets behind them.

It comes as a Python package with three ways in: a **GUI** (Streamlit), a
**command-line tool** and an **importable module**.

## Installation

Requires Python 3.9 or newer.

### From GitHub, straight into an environment

```bash
pip install "ssc-data-finder[gui] @ git+https://github.com/giuso/ssc_data_finder.git"
```

Leave out `[gui]` if you only need the module and the command line.

### From a clone (to run or edit the code)

```bash
git clone https://github.com/giuso/ssc_data_finder.git
cd ssc_data_finder
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

`requirements.txt` lists the same dependencies as `pyproject.toml`, including
the GUI ones. For the notebook, also run `pip install -e ".[notebook]"`.

## GUI

```bash
ssc-data-finder-gui            # opens http://localhost:8501
```

- **Search** (sidebar): a topic, an author, or both. Topics are matched as an
  exact phrase in title/abstract/keywords; write your own query with quotes or
  `AND`/`OR`/`NOT`. Authors can be "John Varga" or "Varga J"; separate several
  with `;`.
- **Filters** (one row above the tables): organism (Human, Mouse, ...),
  publication year, and "only papers with a linked GEO dataset". They apply to
  every table, chart and network on the page.
- **Papers** tab: title, journal, authors, organisms, research centers, GEO
  datasets (clickable), PubMed/DOI links, CSV download.
- **GEO datasets** tab: every GEO series found, including those whose paper did
  not come up in Europe PMC, each with a **Downloads** page listing every file
  (format, size, how to load it in R/Python, `wget` commands, raw reads in SRA).
- **GEO dataset explorer** tab: reads the full GEO record of each dataset (cached
  in `~/.cache/ssc_data_finder/geo` for a week). Charts of samples by tissue,
  condition and sex, datasets by technology, instrument, file format, lab and
  year; a harmonised table of all samples; and a detailed view of any one dataset.
  Samples can be filtered by tissue, condition and technology.
- **Research centers** and **Authors** networks: nodes are institutions/authors,
  linked when they share a paper; size and shade show the number of papers. The
  ranking table next to each network lists papers, collaborators and centrality.
  Author names link to a new search for that author, plus PubMed and Google
  Scholar. A searched author and their centers are highlighted in orange.

The address bar always reflects the current search, so searches can be
bookmarked: `http://localhost:8501/?topic=systemic+sclerosis&author=John+Varga`.
`http://localhost:8501/?downloads=GSE264508` opens the downloads page of any GEO series.

## Command line

```bash
ssc-data-finder --topic "systemic sclerosis"
ssc-data-finder --topic "systemic sclerosis" --max-results 50 --format csv --output results.csv
ssc-data-finder --topic "systemic sclerosis" --format json
ssc-data-finder --author "John Varga"                              # all their scRNA-seq papers
ssc-data-finder --topic "systemic sclerosis" --author "Varga J"    # their papers on the topic
```

## Python module

```python
import pandas as pd
from ssc_data_finder import search_all, run_search, fetch_details, institution_graph, author_graph, ranking

result = search_all("systemic sclerosis", max_results=50)        # also: author="John Varga"
papers = pd.DataFrame(result["papers"])                          # one row per paper
datasets = pd.DataFrame(result["datasets"])                      # one row per GEO series

papers[papers["has_known_dataset"]][["title", "linked_geo_accessions", "organisms"]]

# Most active research centers and authors
pd.DataFrame(ranking(institution_graph(result["papers"]))).head(10)
pd.DataFrame(ranking(author_graph(result["papers"]))).head(10)

# Full GEO record of one dataset: samples, characteristics, files
gse = fetch_details("GSE264508")
gse.title, gse.institute, gse.technologies, gse.file_formats
gse.values("disease")                                            # Counter of samples per condition
[(f.name, f.format, f.url) for f in gse.supplementary_files]
```

`run_search(...)` returns just the paper rows (same as `search_all(...)["papers"]`).

### What "linked dataset" means

Each paper row has the GEO accessions mentioned in its abstract
(`mentioned_geo_accessions`) and the GEO series linked to it
(`linked_geo_accessions`, `has_known_dataset`). A series counts as linked when

- NCBI links it to the paper's PMID (the dataset was deposited with the paper),
- the GEO record lists the paper's PMID,
- the abstract mentions an accession the GEO search also found, or
- the GEO series has (nearly) the same title as the paper.

A paper's `organisms` come from its MeSH tags and its linked datasets, falling
back to words in the abstract ("patients", "mice", ...).

## How it works

- `europepmc.py` — queries the [Europe PMC REST API](https://europepmc.org/RestfulWebService)
  for papers matching the topic/author AND common scRNA-seq phrasing; parses
  authors, affiliations (normalised to institutions) and MeSH organisms.
- `geo.py` — queries [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html)
  (`esearch`/`esummary`/`elink` on the `gds` database) for GEO series matching
  the same search and for the series linked to each paper.
- `geo_details.py` — fetches and parses full GEO series/sample records (SOFT text)
  and FTP folder listings (file sizes).
- `core.py` — runs the searches and cross-references papers and datasets.
- `network.py` — builds the research-center and author graphs (networkx).
- `organisms.py`, `query.py` — shared helpers (organism names, query building).
- `app.py`, `explorer.py`, `downloads.py`, `gui.py` — the Streamlit GUI.
- `cli.py` — the command-line tool.

Institution names are extracted heuristically from free-text affiliations, and
GEO sample characteristics are free text, so the same center or condition can
occasionally appear under two spellings.

NCBI E-utilities allows a higher rate limit with an API key — set
`NCBI_API_KEY` and wire it through if you hit rate limits on larger searches.
