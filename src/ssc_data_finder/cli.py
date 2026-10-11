"""CLI: find published scRNA-seq papers on a topic and cross-reference GEO datasets."""
from __future__ import annotations

import csv
import json
import sys

import click

from ssc_data_finder.core import run_search


@click.command()
@click.option("--topic", default=None, help='Topic/keyword to search, e.g. "systemic sclerosis"')
@click.option("--author", default=None, help='Author name, e.g. "John Varga" or "Varga J"; several separated by ";"')
@click.option("--max-results", default=25, show_default=True, help="Max papers/datasets to fetch per source.")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "json", "csv"]),
    default="table",
    show_default=True,
)
@click.option("--output", type=click.Path(), default=None, help="Write results to this file instead of stdout.")
def main(topic: str | None, author: str | None, max_results: int, output_format: str, output: str | None) -> None:
    """Find papers using single-cell RNA-seq on a topic and/or by an author, and cross-reference GEO datasets."""
    if not topic and not author:
        raise click.UsageError("Give --topic, --author, or both.")
    label = " by ".join(x for x in (topic and f"'{topic}'", author) if x)
    click.echo(f"Searching Europe PMC and GEO for scRNA-seq results on {label}...", err=True)
    rows = run_search(topic, max_results=max_results, author=author)

    if output_format == "json":
        text = json.dumps(rows, indent=2)
    elif output_format == "csv":
        buf = sys.stdout if output is None else open(output, "w", newline="")
        writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "mentioned_geo_accessions": ";".join(row["mentioned_geo_accessions"]),
                    "linked_geo_accessions": ";".join(row["linked_geo_accessions"]),
                    "authors": ";".join(a["name"] for a in row["authors"]),
                    "institutions": ";".join(row["institutions"]),
                    "organisms": ";".join(row["organisms"]),
                }
            )
        if output is not None:
            buf.close()
            click.echo(f"Wrote {len(rows)} rows to {output}", err=True)
        return
    else:
        lines = [f"Found {len(rows)} papers for {label}:\n"]
        for row in rows:
            flag = "[dataset]" if row["has_known_dataset"] else ""
            lines.append(f"- {row['title']} ({row['pub_year']}, {row['journal']}) {flag}")
            if row["pmid"]:
                lines.append(f"    PMID: {row['pmid']}")
            if row["linked_geo_accessions"]:
                lines.append(f"    GEO: {', '.join(row['linked_geo_accessions'])}")
        text = "\n".join(lines)

    if output:
        with open(output, "w") as f:
            f.write(text)
        click.echo(f"Wrote results to {output}", err=True)
    else:
        click.echo(text)


if __name__ == "__main__":
    main()
