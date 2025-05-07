import pandas as pd

# File paths
combined_file = 'combined_motifs.tsv'  # Before filtering
filtered_file = 'filtered_motifs_1000bp_qvalue.tsv'  # After filtering (q-value < 0.05)

# Read the files into pandas DataFrames
combined_df = pd.read_csv(combined_file, sep='\t', header=None, names=['motif_id', 'tf_name', 'gene_id', 'start', 'end', 'strand', 'score', 'p_value', 'q_value', 'motif_sequence'])
filtered_df = pd.read_csv(filtered_file, sep='\t', header=None, names=['motif_id', 'tf_name', 'gene_id', 'start', 'end', 'strand', 'score', 'p_value', 'q_value', 'motif_sequence'])

# Count motifs per gene in the combined file (before filtering)
starting_motifs = combined_df.groupby('gene_id')['motif_id'].count().reset_index()
starting_motifs.columns = ['gene_id', 'starting_motifs']

# Count motifs per gene in the filtered file (after applying q-value < 0.05 filter)
filtered_motifs = filtered_df.groupby('gene_id')['motif_id'].count().reset_index()
filtered_motifs.columns = ['gene_id', 'filtered_motifs']

# Merge both DataFrames to create a results table
result_table = pd.merge(starting_motifs, filtered_motifs, on='gene_id', how='left').fillna(0)

# Save the results table to a CSV file
result_table.to_csv('motif_counts_per_gene.csv', index=False)

# Display the results
print(result_table)


# Load the filtered motif file
df = pd.read_csv("filtered_motifs_1000bp_qvalue.tsv", sep="\t", header=None)

# Add column names for easier understanding
df.columns = ["motif_id", "motif_name", "gene_id", "start", "end", "strand", "score", "p-value", "q-value",
              "matched_sequence"]


# Function to extract gene group/type (e.g., 'melanocyte', 'SRY-box') from gene_id column
def extract_gene_group(gid):
    if "|" in gid:
        return gid.split("|")[-1].strip()  # Get the last part (e.g. 'melanocyte')
    else:
        return "unknown"



# Add gene group/type column
df["gene_group"] = df["gene_id"].astype(str).apply(extract_gene_group)


# Group by gene group
grouped = df.groupby("gene_group")

# Compare motif IDs between all genes in the same group
shared_motifs = {}

for group_name, group_df in grouped:
    genes = group_df["gene_id"].unique()
    motif_sets = {gene: set(group_df[group_df["gene_id"] == gene]["motif_id"]) for gene in genes}

    if len(motif_sets) > 1:
        # Take intersection of all motif sets for this group
        common_motifs = set.intersection(*motif_sets.values())
        shared_motifs[group_name] = {
            "genes": list(genes),
            "shared_motifs": common_motifs,
            "count": len(common_motifs)
        }

# Print results
for group, data in shared_motifs.items():
    print(f"\nGene Group: {group}")
    print(f"Genes: {data['genes']}")
    print(f"Shared Motifs Count: {data['count']}")
    print(f"Shared Motifs: {sorted(data['shared_motifs'])}")

# Convert shared_motifs dictionary into a DataFrame
output_rows = []
for group, data in shared_motifs.items():
    output_rows.append({
        "Gene Group": group,
        "Orthologous Genes": ", ".join(data["genes"]),
        "Shared Motif Count": data["count"],
        "Shared Motifs": ", ".join(sorted(data["shared_motifs"]))
    })

shared_motifs_df = pd.DataFrame(output_rows)

# Save to CSV
shared_motifs_df.to_csv("shared_conserved_motifs.csv", index=False)

# Optional: Print table nicely (Markdown format)
print("\nConserved Motifs Table:")
print(shared_motifs_df)

