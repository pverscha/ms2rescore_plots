import os
import requests
import pandas as pd
from tqdm import tqdm
from collections import Counter

def query_api(url, payload):
    """
    General function to perform API requests.
    """
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=payload, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"API request failed with status code {response.status_code}: {response.text}")

def get_taxa_info(taxids):
    """
    Retrieves taxa information, including lineage, for a list of tax IDs.
    """
    url = "https://api.unipept.ugent.be/private_api/taxa"
    payload = {"taxids": list(filter(None, taxids))}  # Remove None values
    batch_size = 100
    all_taxa = {}
    for i in range(0, len(payload["taxids"]), batch_size):
        batch = payload["taxids"][i:i + batch_size]
        response = query_api(url, {"taxids": batch})
        for taxon in response:
            all_taxa[taxon["id"]] = taxon
    return all_taxa

def process_files(input_directory, output_directory):
    os.makedirs(output_directory, exist_ok=True)
    
    for filename in os.listdir(input_directory):
        if filename.endswith(".txt"):
            input_filepath = os.path.join(input_directory, filename)
            output_filepath = os.path.join(output_directory, filename.replace(".txt", "_taxonomy.tsv"))
            
            with open(input_filepath, 'r') as file:
                peptides = [line.strip() for line in file.readlines()]

            peptide_responses = []
            all_taxids = set()
            
            # Process peptides in batches of 10
            with tqdm(total=len(peptides), desc=f"Processing {filename}", unit="peptide") as pbar:
                for i in range(0, len(peptides), 10):
                    batch = peptides[i:i+10]
                    response_data = query_api("https://api.unipept.ugent.be/mpa/pept2data", {"peptides": batch, "equate_il": True, "missed": True})
                    for item in response_data.get("peptides", []):
                        lca = item.get("lca")
                        if lca:
                            all_taxids.add(lca)
                            peptide_responses.append((item["sequence"], lca))
                    pbar.update(len(batch))

            # Retrieve the LCA taxonomies
            taxid_to_info = get_taxa_info(all_taxids)
            # Collect all lineage taxids for a further lookup
            all_lineage_taxids = set()
            for taxid, info in taxid_to_info.items():
                if info.get("lineage"):
                    all_lineage_taxids.update(info["lineage"])

            # Retrieve detailed lineage names
            detailed_taxid_to_name = get_taxa_info(all_lineage_taxids)

            # Prepare final data
            final_data = []
            for peptide, lca in peptide_responses:
                taxa_info = taxid_to_info.get(lca, {})
                lineage = taxa_info.get("lineage", [None]*29)
                lineage_names = [detailed_taxid_to_name.get(taxid, {}).get("name", "") for taxid in lineage]
                final_data.append((peptide, taxa_info.get("name", "")) + tuple(lineage_names))

            # Define the columns
            columns = ["peptide", "lca", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "superorder", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species group", "species subgroup", "species", "subspecies", "strain", "varietas", "forma"]
            df = pd.DataFrame(final_data, columns=columns)
            df.to_csv(output_filepath, sep='\t', index=False)
            print(f"Saved results to {output_filepath}")

# Set the input and output directories
input_directory = './peptides'
output_directory = './lcas'

# Run the processing function
process_files(input_directory, output_directory)
