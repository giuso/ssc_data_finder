"""Find published single-cell RNA-seq papers and their GEO datasets.

Main entry points:

    from ssc_data_finder import search_all, run_search, fetch_details

    result = search_all("systemic sclerosis", max_results=50)   # {"papers": [...], "datasets": [...]}
    papers = run_search(author="John Varga")                      # paper rows only
    details = fetch_details("GSE264508")                          # full GEO record + samples
"""
from ssc_data_finder.core import run_search, search_all
from ssc_data_finder.geo import geo_url
from ssc_data_finder.geo_details import SeriesDetails, fetch_details, fetch_many
from ssc_data_finder.network import author_graph, institution_graph, ranking

__version__ = "0.1.0"

__all__ = [
    "search_all",
    "run_search",
    "fetch_details",
    "fetch_many",
    "SeriesDetails",
    "geo_url",
    "institution_graph",
    "author_graph",
    "ranking",
]
