"""Streamlit GUI: results table, then research-center and author collaboration networks.

Run with `ssc-data-finder-gui` (or `streamlit run src/ssc_data_finder/app.py`).
"""
from __future__ import annotations

from urllib.parse import quote, urlencode

import networkx as nx
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from pyvis.network import Network

from ssc_data_finder import downloads, explorer, network
from ssc_data_finder.organisms import NOT_STATED, ordered
from ssc_data_finder.core import search_all
from ssc_data_finder.geo import geo_url
from ssc_data_finder.query import author_matches, last_initials, split_authors

# Sequential blue ramp (light -> dark), steps 250..650 of the reference palette.
# Darker = more papers; the lightest step still clears 2:1 on the light surface.
NODE_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
NODE_RAMP_DARK = ["#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4"]
SURFACE = {"light": "#fcfcfb", "dark": "#1a1a19"}
TEXT = {"light": "#0b0b0b", "dark": "#f0efec"}
EDGE = {"light": "#c3c2b7", "dark": "#52514e"}
# Categorical slot 2 (orange): the searched author(s) and their research centers.
HIGHLIGHT = {"light": "#eb6834", "dark": "#d95926"}
LABELLED_NODES = 12

st.set_page_config(page_title="scRNA-seq Data Finder", layout="wide")


@st.cache_data(ttl=3600, show_spinner=False)
def search(topic: str, author: str, max_results: int) -> dict:
    return search_all(topic or None, max_results=max_results, author=author or None)


def geo_link(accessions: list[str]) -> str | None:
    """One clickable cell for any number of GEO accessions.

    A single accession links to its GEO page, several to a GEO search listing them all.
    The "#GSE1, GSE2" fragment is ignored by NCBI and only feeds the cell's display text.
    """
    if not accessions:
        return None
    if len(accessions) == 1:
        url = geo_url(accessions[0])
    else:
        url = "https://www.ncbi.nlm.nih.gov/gds/?term=" + quote(" OR ".join(f"{a}[ACCN]" for a in accessions))
    return f"{url}#{', '.join(accessions)}"


def pubmed_link(pmids: list[str]) -> str | None:
    if not pmids:
        return None
    if len(pmids) == 1:
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/"
    else:
        url = "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote(" OR ".join(f"{p}[uid]" for p in pmids))
    return f"{url}#{', '.join(pmids)}"


LINK_TEXT = r"#(.+)$"  # display the fragment written by geo_link/pubmed_link


def theme() -> str:
    try:
        return st.context.theme.type or "light"
    except AttributeError:
        return "light"


def results_table(rows: list[dict]) -> pd.DataFrame:
    def short_authors(authors: list[dict]) -> str:
        names = [a["name"] for a in authors]
        return ", ".join(names[:3]) + (" et al." if len(names) > 3 else "")

    return pd.DataFrame(
        {
            "Title": [r["title"] for r in rows],
            "Year": [r["pub_year"] for r in rows],
            "Journal": [r["journal"] or "" for r in rows],
            "Authors": [short_authors(r["authors"]) for r in rows],
            "Organism": [", ".join(r["organisms"]) for r in rows],
            "Research centers": ["; ".join(r["institutions"]) for r in rows],
            "GEO datasets": [geo_link(r["linked_geo_accessions"]) for r in rows],
            "GEO in abstract": [geo_link(r["mentioned_geo_accessions"]) for r in rows],
            "Dataset": [r["has_known_dataset"] for r in rows],
            "PubMed": [f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/" if r["pmid"] else None for r in rows],
            "DOI": [f"https://doi.org/{r['doi']}" if r["doi"] else None for r in rows],
        }
    )


def organism_filter(papers: list[dict], datasets: list[dict]) -> set[str] | None:
    """Organism pills, most common first. Returns the selected set, or None when nothing is excluded."""
    counts: dict[str, int] = {}
    for item in papers + datasets:
        for organism in item["organisms"] or [NOT_STATED]:
            counts[organism] = counts.get(organism, 0) + 1
    options = ordered(o for o in counts if o != NOT_STATED) + ([NOT_STATED] if NOT_STATED in counts else [])
    selected = st.pills(
        "Organism",
        options,
        selection_mode="multi",
        default=options,
        key="organisms_" + "_".join(options),  # reset when a new search brings different organisms
        help="Papers and datasets are kept if they involve any selected organism. "
        "Paper organisms come from MeSH tags, linked GEO datasets, or the abstract.",
    )
    if not selected or set(selected) == set(options):
        return None
    return set(selected)


def searched_nodes(rows: list[dict], author: str) -> tuple[set[str], set[str]]:
    """Author names matching the searched author(s), and every research center they list."""
    queries = split_authors(author)
    names, centers = set(), set()
    for row in rows:
        for a in row["authors"]:
            if any(author_matches(a["name"], q) for q in queries):
                names.add(a["name"])
                centers.update(a["institutions"])
    return names, centers


def render_graph(graph: nx.Graph, mode: str, highlight: set[str], height: int = 620) -> None:
    ramp = NODE_RAMP_DARK if mode == "dark" else NODE_RAMP
    max_papers = max(d["papers"] for _, d in graph.nodes(data=True))
    labelled = {r["name"] for r in network.ranking(graph)[:LABELLED_NODES]} | highlight

    net = Network(
        height=f"{height}px",
        width="100%",
        bgcolor=SURFACE[mode],
        font_color=TEXT[mode],
        cdn_resources="remote",
    )
    for node, data in graph.nodes(data=True):
        step = round((data["papers"] - 1) / max(max_papers - 1, 1) * (len(ramp) - 1))
        tooltip = [node, f"{data['papers']} paper(s) · {graph.degree(node)} collaborator(s)"]
        if data.get("institution"):
            tooltip.append(data["institution"])
        tooltip += [f"• {t[:90]}" for t in data["titles"][:5]]
        if len(data["titles"]) > 5:
            tooltip.append(f"… and {len(data['titles']) - 5} more")
        net.add_node(
            node,
            label=node if node in labelled else " ",
            title="\n".join(tooltip),
            size=8 + 22 * (data["papers"] / max_papers) ** 0.5,
            color={
                "background": HIGHLIGHT[mode] if node in highlight else ramp[step],
                "border": SURFACE[mode],
                "highlight": HIGHLIGHT[mode] if node in highlight else ramp[-1],
            },
            borderWidth=2,
        )
    for a, b, data in graph.edges(data=True):
        net.add_edge(a, b, value=data["weight"], color=EDGE[mode], title=f"{data['weight']} shared paper(s)")
    net.barnes_hut(gravity=-6000, spring_length=140, overlap=0.4)
    net.set_edge_smooth("continuous")
    components.html(net.generate_html(), height=height + 20)


def app_search_url(author: str, topic: str = "") -> str:
    """Link back into this app that runs a new search (read from st.query_params on load)."""
    return "?" + urlencode({"author": author, "topic": topic}) + f"#{author}"


def author_link_columns(ranking: pd.DataFrame, topic: str) -> pd.DataFrame:
    names = list(ranking["name"])
    links = {
        "name": [app_search_url(n) for n in names],
        "in_topic": [app_search_url(n, topic) for n in names] if topic else None,
        "pubmed": ["https://pubmed.ncbi.nlm.nih.gov/?term=" + quote(f'"{last_initials(n)}"[Author]') for n in names],
        "scholar": ["https://scholar.google.com/scholar?q=" + quote(f'author:"{n}"') for n in names],
    }
    return ranking.assign(**{k: v for k, v in links.items() if v is not None})


def network_section(
    title: str,
    caption: str,
    graph: nx.Graph,
    default_n: int,
    key: str,
    mode: str,
    highlight: set[str],
    author_links: bool = False,
) -> None:
    st.subheader(title)
    st.caption(caption)
    if graph.number_of_nodes() == 0:
        st.info("No affiliation data for these papers.")
        return
    n = st.slider(
        "Nodes shown (most papers first)",
        min_value=min(5, graph.number_of_nodes()),
        max_value=graph.number_of_nodes(),
        value=min(default_n, graph.number_of_nodes()),
        key=key,
    )
    # Searched author(s) and their centers are always shown, whatever the slider says.
    shown = graph.subgraph(set(network.top_nodes(graph, n)) | (highlight & set(graph))).copy()
    ranking = pd.DataFrame(network.ranking(graph))
    if highlight:
        ranking.insert(0, "searched", ranking["name"].isin(highlight))
    if author_links:
        ranking = author_link_columns(ranking, st.session_state.topic)
    order = [c for c in ["searched", "name", "in_topic", "pubmed", "scholar", "papers", "collaborators",
                         "centrality", "institution"] if c in ranking]
    left, right = st.columns([2, 1])
    with left:
        render_graph(shown, mode, highlight)
    with right:
        st.dataframe(
            ranking,
            hide_index=True,
            height=620,
            column_order=order,
            column_config={
                "searched": st.column_config.CheckboxColumn("★", help="Searched author / their research center"),
                "name": st.column_config.LinkColumn(
                    "Name", display_text=r"#(.+)$", help="Opens a new search for this author (all topics)"
                )
                if author_links
                else "Name",
                "in_topic": st.column_config.LinkColumn(
                    "Search", display_text="in topic", help="This author's papers on the current topic"
                ),
                "pubmed": st.column_config.LinkColumn("PubMed", display_text="PubMed"),
                "scholar": st.column_config.LinkColumn("Scholar", display_text="Scholar"),
                "papers": st.column_config.NumberColumn("Papers"),
                "collaborators": st.column_config.NumberColumn("Collaborators", help="Distinct co-author nodes"),
                "centrality": st.column_config.ProgressColumn(
                    "Centrality", help="Share of the network directly connected", min_value=0.0, max_value=1.0
                ),
                "institution": "Main research center",
            },
        )


# --- Downloads page: "?downloads=GSExxxx" (opened from the GEO tables) ---------
if st.query_params.get("downloads"):
    downloads.render(st.query_params["downloads"])
    st.stop()

# --- Sidebar: search controls ------------------------------------------------
with st.sidebar:
    st.header("Search")
    with st.form("search"):
        topic = st.text_input("Topic", value=st.query_params.get("topic", "systemic sclerosis"), help='Exact phrase; use quotes or AND/OR/NOT for your own query')
        author = st.text_input(
            "Author",
            value=st.query_params.get("author", ""),
            placeholder="e.g. John Varga or Varga J",
            help='Optional. Several authors separated by ";". With a topic, finds that author\'s papers on the topic.',
        )
        max_results = st.number_input("Max results per source", min_value=5, max_value=500, value=50, step=5)
        submitted = st.form_submit_button("Search", type="primary", width="stretch")

if submitted and not topic.strip() and not author.strip():
    st.sidebar.error("Enter a topic, an author, or both.")
    submitted = False

if submitted or "rows" not in st.session_state:
    label = " by ".join(x for x in (topic.strip() and f"'{topic.strip()}'", author.strip()) if x)
    with st.spinner(f"Searching Europe PMC and GEO for {label}…"):
        try:
            result = search(topic.strip(), author.strip(), int(max_results))
            st.session_state.rows = result["papers"]
            st.session_state.datasets = result["datasets"]
            st.session_state.topic = topic.strip()
            st.session_state.author = author.strip()
            # Keep the address bar in sync so a search can be bookmarked or shared.
            st.query_params.from_dict({"topic": topic.strip(), "author": author.strip()})
        except Exception as exc:  # network/API failures shouldn't kill the page
            st.error(f"Search failed: {exc}")
            st.stop()

all_rows = st.session_state.rows
all_datasets = st.session_state.datasets
mode = theme()

st.title("scRNA-seq Data Finder")
searched = [f"Topic: **{st.session_state.topic}**"] if st.session_state.topic else []
if st.session_state.author:
    searched.append(f"Author: **{st.session_state.author}**")
st.caption(" · ".join(searched + ["Europe PMC papers cross-referenced with NCBI GEO series"]))

# --- Filters: one row, applied to every table, chart and network below ---------
f_org, f_year, f_geo = st.columns([3, 3, 2], vertical_alignment="bottom")
with f_org:
    selected_organisms = organism_filter(all_rows, all_datasets)
with f_year:
    years = sorted({int(r["pub_year"]) for r in all_rows if str(r["pub_year"] or "").isdigit()})
    if len(years) > 1:
        year_range = st.slider("Publication year", years[0], years[-1], (years[0], years[-1]), key=f"years_{years[0]}_{years[-1]}")
    else:
        year_range = None
with f_geo:
    only_datasets = st.toggle("Only papers with a linked GEO dataset", value=False)


def keep_paper(r: dict) -> bool:
    if selected_organisms is not None and not set(r["organisms"] or [NOT_STATED]) & selected_organisms:
        return False
    if year_range and str(r["pub_year"] or "").isdigit() and not year_range[0] <= int(r["pub_year"]) <= year_range[1]:
        return False
    return r["has_known_dataset"] or not only_datasets


rows = [r for r in all_rows if keep_paper(r)]
datasets = [
    d for d in all_datasets if selected_organisms is None or set(d["organisms"] or [NOT_STATED]) & selected_organisms
]

centers = network.institution_graph(rows)
authors = network.author_graph(rows)
searched_authors, searched_centers = searched_nodes(rows, st.session_state.author)
highlight_note = " Orange = the searched author(s) and their research centers." if searched_authors else ""

c1, c2, c3, c4 = st.columns(4)
c1.metric("Papers", len(rows))
c2.metric("With linked GEO dataset", sum(r["has_known_dataset"] for r in rows))
c3.metric("Research centers", centers.number_of_nodes())
c4.metric("Authors", authors.number_of_nodes())

# --- Results tables ----------------------------------------------------------
file_stem = "_".join(filter(None, [st.session_state.topic, st.session_state.author])).replace(" ", "_")
papers_tab, datasets_tab, explorer_tab = st.tabs(
    [f"Papers ({len(rows)})", f"GEO datasets ({len(datasets)})", "GEO dataset explorer"]
)

with papers_tab:
    df = results_table(rows)
    st.dataframe(
        df,
        hide_index=True,
        width="stretch",
        column_config={
            "Title": st.column_config.TextColumn(width="large"),
            "Year": st.column_config.TextColumn(width="small"),
            "GEO datasets": st.column_config.LinkColumn(
                display_text=LINK_TEXT, help="GEO series deposited with this paper (linked by NCBI or confirmed by the GEO search)"
            ),
            "GEO in abstract": st.column_config.LinkColumn(
                display_text=LINK_TEXT, help="Accessions written in the abstract (may include re-used public data)"
            ),
            "Dataset": st.column_config.CheckboxColumn(help="At least one GEO series is linked to this paper"),
            "PubMed": st.column_config.LinkColumn(display_text="PubMed"),
            "DOI": st.column_config.LinkColumn(display_text="DOI"),
        },
    )
    export = df.assign(
        **{
            "GEO datasets": [", ".join(r["linked_geo_accessions"]) for r in rows],
            "GEO in abstract": [", ".join(r["mentioned_geo_accessions"]) for r in rows],
        }
    )
    st.download_button(
        "Download CSV", export.to_csv(index=False).encode(), file_name=f"{file_stem}_scrnaseq_papers.csv", mime="text/csv"
    )

with datasets_tab:
    st.caption(
        "Every GEO series found by the GEO search or linked by NCBI to one of the papers. "
        "Series whose paper is not in the Papers tab are papers found through GEO."
    )
    ds = pd.DataFrame(
        {
            "Accession": [d["url"] for d in datasets],
            "Title": [d["title"] for d in datasets],
            "Organism": [", ".join(d["organisms"]) for d in datasets],
            "Samples": [d["n_samples"] for d in datasets],
            "Paper (PubMed)": [pubmed_link(d["pubmed_ids"]) for d in datasets],
            "Paper in Papers tab": [d["paper_in_results"] for d in datasets],
            "Downloads": [downloads.downloads_url(d["accession"]) for d in datasets],
        }
    )
    st.dataframe(
        ds,
        hide_index=True,
        width="stretch",
        column_config={
            "Accession": st.column_config.LinkColumn(display_text=r"acc=(GSE\d+)"),
            "Title": st.column_config.TextColumn(width="large"),
            "Paper (PubMed)": st.column_config.LinkColumn(
                display_text=LINK_TEXT, help="Papers the GEO record cites; empty if none published yet"
            ),
            "Paper in Papers tab": st.column_config.CheckboxColumn(),
            "Downloads": st.column_config.LinkColumn(
                display_text="files", help="Opens a page listing every file you can download, with sizes"
            ),
        },
    )
    st.download_button(
        "Download CSV",
        ds.assign(
            Accession=[d["accession"] for d in datasets],
            Downloads=[d["url"].replace("query/acc.cgi", "download/") for d in datasets],
            **{"Paper (PubMed)": [", ".join(d["pubmed_ids"]) for d in datasets]},
        ).to_csv(index=False).encode(),
        file_name=f"{file_stem}_geo_datasets.csv",
        mime="text/csv",
        key="datasets_csv",
    )

with explorer_tab:
    explorer.render(datasets, mode, selected_organisms)

# --- Networks ----------------------------------------------------------------
network_section(
    "Research centers",
    "Each node is an institution; an edge means the two co-authored a paper. "
    "Bigger and darker = more papers. Hover for details, drag to rearrange, scroll to zoom." + highlight_note,
    centers,
    default_n=40,
    key="centers_n",
    mode=mode,
    highlight=searched_centers,
)
network_section(
    "Authors",
    "Each node is an author; an edge means they share a paper. "
    "Bigger and darker = more papers. Hover for the author's main research center and papers. "
    "In the table, click a name to search that author (opens a new tab)." + highlight_note,
    authors,
    default_n=60,
    key="authors_n",
    mode=mode,
    highlight=searched_authors,
    author_links=True,
)
