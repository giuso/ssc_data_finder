"""Build collaboration networks (research centers, authors) from search result rows."""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations

import networkx as nx


def institution_graph(rows: list[dict]) -> nx.Graph:
    """Institutions as nodes; an edge joins two institutions that co-authored a paper."""
    graph = nx.Graph()
    papers = defaultdict(list)
    for row in rows:
        for inst in row["institutions"]:
            papers[inst].append(row["title"])
        for a, b in combinations(row["institutions"], 2):
            _add_edge(graph, a, b)
    for inst, titles in papers.items():
        graph.add_node(inst, papers=len(titles), titles=titles)
    return graph


def author_graph(rows: list[dict]) -> nx.Graph:
    """Authors as nodes; an edge joins two authors who share a paper."""
    graph = nx.Graph()
    papers = defaultdict(list)
    institutions = defaultdict(Counter)
    for row in rows:
        names = list(dict.fromkeys(a["name"] for a in row["authors"]))
        for author in row["authors"]:
            institutions[author["name"]].update(author["institutions"])
        for name in names:
            papers[name].append(row["title"])
        for a, b in combinations(names, 2):
            _add_edge(graph, a, b)
    for name, titles in papers.items():
        top = institutions[name].most_common(1)
        graph.add_node(
            name,
            papers=len(titles),
            titles=titles,
            institution=top[0][0] if top else None,
        )
    return graph


def top_nodes(graph: nx.Graph, n: int) -> nx.Graph:
    """Keep the `n` nodes with most papers (ties broken by number of collaborators)."""
    ranked = sorted(graph.nodes, key=lambda v: (graph.nodes[v]["papers"], graph.degree(v)), reverse=True)
    return graph.subgraph(ranked[:n]).copy()


def ranking(graph: nx.Graph) -> list[dict]:
    """One row per node: paper count, collaborator count and centrality, best first."""
    centrality = nx.degree_centrality(graph) if len(graph) > 1 else {v: 0.0 for v in graph}
    out = [
        {
            "name": v,
            "papers": data["papers"],
            "collaborators": graph.degree(v),
            "centrality": round(centrality[v], 3),
            **({"institution": data.get("institution")} if "institution" in data else {}),
        }
        for v, data in graph.nodes(data=True)
    ]
    return sorted(out, key=lambda r: (r["papers"], r["collaborators"]), reverse=True)


def _add_edge(graph: nx.Graph, a: str, b: str) -> None:
    if graph.has_edge(a, b):
        graph[a][b]["weight"] += 1
    else:
        graph.add_edge(a, b, weight=1)
