#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import glob
import json
import os
import sys
from typing import Any

from eutils import Client
from mibig.converters.shared.common import Citation

from mibig_html.annotations.references import PubmedCache, PubmedEntry


def gather_references(data: Any, existing_results: list[str] = None) -> list[Citation]:
    if existing_results is None:
        existing_results = []
    if isinstance(data, dict):
        for key, val in data.items():
            gather_references(val, existing_results)
    elif isinstance(data, str):
        if data.startswith("pubmed:"):
            data = data.split()[0]
            existing_results.append(Citation.from_json(data))
    elif isinstance(data, list):
        for val in data:
            gather_references(val, existing_results)
    return existing_results


def extract_pmids(files: list[str]) -> list[str]:
    """Extract all unique pmids from mibig json files"""
    for filename in files:
        with open(filename) as handle:
            data = json.load(handle)
        citations = gather_references(data)
    return sorted(list({citation.value for citation in citations}))


def fetch_all(cache_file: str, citations: list[Citation]) -> None:
    pubmed_cache = PubmedCache(cache_file)
    client = Client(api_key=os.environ.get("NCBI_API_KEY", None))

    missing = pubmed_cache.get_missing(pmids)
    while missing:
        articles = client.efetch(db="pubmed", id=missing)
        for article in articles:
            pubmed_cache.add(article.title, article.authors,
                                article.year, article.jrnl, article.pmid)
        missing = pubmed_cache.get_missing(pmids)
        
    pubmed_cache.save()


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} input_dir cache_file")
        sys.exit(1)
    cache = "pubmed_cache.json"
    if len(sys.argv) == 3:
        cache = sys.argv[2]
    pmids = extract_pmids(glob.glob(f"{sys.argv[1]}/*.json"))
    fetch_all(cache, pmids)
