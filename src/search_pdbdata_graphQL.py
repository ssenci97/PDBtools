#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import requests
import pandas as pd

def search_rcsb(query_text, min_year):
    """Query RCSB Search API v2 combining a free-text search with a release year limit."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    all_pdb_ids = []
    rows_per_page = 10000
    start = 0

    # Build a grouped query: Full-text match AND Released >= Jan 1st of min_year
    search_payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": query_text
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "greater_or_equal",
                        "value": f"{min_year}-01-01"
                    }
                }
            ]
        },
        "request_options": {
            "paginate": {
                "start": start,
                "rows": rows_per_page
            }
        },
        "return_type": "entry"
    }

    while True:
        try:
            response = requests.post(
                url,
                json=search_payload,
                headers={"Content-Type": "application/json"},
                timeout=60
            )
            response.raise_for_status()
            if not response.text:
                break
            data = response.json()

            total_count = data.get("total_count", 0)
            result_set = data.get("result_set", [])
            page_ids = [entry["identifier"] for entry in result_set]
            all_pdb_ids.extend(page_ids)

            if len(all_pdb_ids) >= total_count or len(page_ids) == 0:
                break

            start += rows_per_page
            search_payload["request_options"]["paginate"]["start"] = start

        except Exception as e:
            print(f"Warning: Search failed at offset {start}: {e}")
            break

    return all_pdb_ids


def fetch_metadata(pdb_ids):
    """Fetch entry-level metadata using GraphQL for the matching PDB IDs."""
    if not pdb_ids:
        return pd.DataFrame()
        
    url = "https://data.rcsb.org/graphql"
    
    graphql_query = """
    query GetPdbDetails($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        rcsb_accession_info {
          initial_release_date
        }
        rcsb_entry_info {
          resolution_combined
        }
        exptl {
          method
        }
        struct {
          title
        }
      }
    }
    """
    
    # We batch in groups of 100 to prevent hitting API payload limits on giant search result sets
    batch_size = 100
    all_rows = []
    
    for i in range(0, len(pdb_ids), batch_size):
        batch = pdb_ids[i:i + batch_size]
        payload = {
            "query": graphql_query,
            "variables": {"ids": batch}
        }
        try:
            response = requests.post(url, json=payload, timeout=60)
            response.raise_for_status()
            entries = response.json().get("data", {}).get("entries") or []
            
            for entry in entries:
                pdb_id = entry.get("rcsb_id")
                release_date = entry.get("rcsb_accession_info", {}).get("initial_release_date")
                resolutions = entry.get("rcsb_entry_info", {}).get("resolution_combined")
                resolution = resolutions[0] if resolutions else None
                exptl = entry.get("exptl") or []
                method = exptl[0].get("method") if exptl else None
                title = entry.get("struct", {}).get("title")
                
                all_rows.append({
                    "pdb_id": pdb_id,
                    "release_date": release_date,
                    "resolution": resolution,
                    "method": method,
                    "title": title
                })
        except Exception as e:
            print(f"Warning: Metadata fetch failed for batch {i}: {e}")
            
    return pd.DataFrame(all_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Search RCSB with full-text queries and fetch entry metadata without downloading structures."
    )
    parser.add_argument("--query", required=True, help="Your search query (e.g. 'kinase', 'ribosome', 'P00533')")
    parser.add_argument("--min_year", type=int, default=2023, help="Keep PDBs released from this year onward (default: 2023)")
    parser.add_argument("--outdir", default="data/pdbs_query", help="Output directory")
    args = parser.parse_args()

    # Sanitize query name for filesystem-safe output filenames
    safe_query_name = "".join(c if c.isalnum() else "_" for c in args.query).strip("_")
    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, f"{safe_query_name}_pdbs.tsv")

    # 1. Get PDB IDs
    print(f"Searching RCSB for: '{args.query}' (Released {args.min_year}+)...")
    pdb_ids = search_rcsb(args.query, args.min_year)
    if not pdb_ids:
        print("No matching structures found.")
        return
    print(f"Found {len(pdb_ids)} structures.")

    # 2. Fetch Metadata
    print("Fetching metadata details via GraphQL...")
    df = fetch_metadata(pdb_ids)

    # Save output
    if not df.empty:
        df.to_csv(out_file, sep="\t", index=False)
        print(f"Saved PDB metadata to {out_file}")
        
        # Display preview
        print(f"\n--- Preview of Results ({len(df)} total) ---")
        print(df[["pdb_id", "release_date", "resolution", "method"]].head(10))
    else:
        print("Failed to compile metadata.")

if __name__ == "__main__":
    main()
