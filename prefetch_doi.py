#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import glob
from json import load
import os
import sys
from typing import Any, List

from mibig.converters.shared.common import Citation

from mibig_html.annotations.references import DoiCache, DoiEntry

SPECIAL = {
    "10.12211/2096-8280.2021-024": {  # times out for anything but HTML
        "title": "Genome mining for novel natural products in Sorangium cellulosum So0157-2 by heterologous expression",
        "authors": ["Zhou, H" "Shen, Q", "Chen, H", "Wang, Z", "Li, Y", "Zhang, Y", "Bian, X"],
        "year": "2021",
        "journal": "Synthetic Biology Journal",
        "identifier": "10.12211/2096-8280.2021-024",
    },
}

def gather_dois(data: Any, existing_results: list[str] = None) -> list[Citation]:
    if existing_results is None:
        existing_results = []
    if isinstance(data, dict):
        for key, val in data.items():
            gather_dois(val, existing_results)
    elif isinstance(data, str):
        if data.startswith("doi:"):
            existing_results.append(Citation.from_json(data))
    elif isinstance(data, list):
        for val in data:
            gather_dois(val, existing_results)
    return existing_results
            

def fetch_all(cache_file: str, files: List[str]) -> None:
    doi_cache = DoiCache(cache_file)
    for filename in files:
        with open(filename) as handle:
            data = load(handle)
        citations = gather_dois(data)
        for citation in citations:
            if citation.value in SPECIAL:
                doi_cache.add_entry(DoiEntry.from_json(SPECIAL[citation.value]))
            else:
                try:
                    doi_cache.get(citation.value)
                    doi_cache.save()
                except ValueError as err:
                    print("failed to import DOIs from", filename, err)
                    raise
    doi_cache.save()


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} input_dir [cache_file]")
        sys.exit(1)
    cache = "doi_cache.json"
    if len(sys.argv) == 3:
        cache = sys.argv[2]
    fetch_all(cache, glob.glob(f"{sys.argv[1]}/*.json"))
