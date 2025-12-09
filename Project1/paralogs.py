#!/usr/bin/env python3

import csv
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

faa_file_path = "/farmshare/home/classes/bios/270/data/project1/ecoli_bakta_out/assembly.faa"
cluster_file_path = "/farmshare/home/classes/bios/270/data/project1/ecoli_mmseqs_out/ecoli_prot90_cluster.tsv"
tsv_out = "paralogs.tsv"
png_out = "paralogs_top10.png"

# functions

def protein_names(faa):
    protein_info = {} # make dictionary

    with open(faa, "r") as names:
        for line in names:
            if line.startswith(">"): # new seq line in fasta bvegins with >
                header = line[1:].strip()
                parts = header.split()
                protein_id = parts[0] # first part of fasta
                
                if len(parts) > 1:
                    protein_name = " ".join(parts[1:]) # second part
                
                else:
                    protein_name = "None" # add value none iof theres no name
                
                protein_info[protein_id] = protein_name

    return protein_info

def load_clusters(cluster_file):
    clusters = defaultdict(list)

    with open(cluster_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")

        for row in reader:
            cluster_id, protein_id = row[0], row[1]
            clusters[cluster_id].append(protein_id)

    return clusters

def build_paralog_table(protein_info, clusters):
    
    rows = []
    
    for cluster_id, members in clusters.items():
        if len(members) <= 1: # only keep groups with multiple paralogs
            continue

        copy_number = len(members)
        for pid in members:
            pname = protein_info.get(pid, "")
            rows.append({
                "protein_id": pid,
                "protein_name": pname,
                "copy_number": copy_number
            })

    if len(rows) == 0:

        print("no paralogs")

    return pd.DataFrame(rows)

def plot_top_paralogs(df, out_png, top_n=10):

    df_unique = df.drop_duplicates(subset=["proteinid"])
    top = df_unique.sort_values("copy_number", ascending=False).head(top_n)

    labels = [
        row["protein_name"] if row["protein_name"] else row["protein_id"]
        for _, row in top.iterrows()
    ]

    plt.figure(figsize=(10, 5))
    plt.bar(range(len(top)), top["copynumber"])
    plt.xticks(range(len(top)), labels, rotation=45, ha="right")
    plt.ylabel("copy number")
    plt.title("paralogs")
    plt.tight_layout()
    plt.savefig(out_png, dpi=600)
    plt.close()

# plotting and running

protein_info = protein_names(faa_file_path)
clusters = load_clusters(cluster_file_path)
df = build_paralog_table(protein_info, clusters)
df.to_csv(tsv_out, sep="\t", index=False)
plot_top_paralogs(df, png_out)
