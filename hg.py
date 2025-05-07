import pandas as pd

# Load data
combined_df = pd.read_csv('combined_motifs.tsv', sep='\t', header=None,
                          names=['motif_id', 'tf_name', 'gene_id', 'start', 'end', 'strand', 'score', 'p_value', 'q_value', 'motif_sequence'])
filtered_df = pd.read_csv('filtered_motifs_1000bp_qvalue.tsv', sep='\t', header=None,
                          names=['motif_id', 'tf_name', 'gene_id', 'start', 'end', 'strand', 'score', 'p_value', 'q_value', 'motif_sequence'])

# Extract gene group from gene_id
def extract_gene_group(gid):
    if "|" in gid:
        return gid.split("|")[-1].strip()
    return "unknown"

combined_df['gene_group'] = combined_df['gene_id'].astype(str).apply(extract_gene_group)
filtered_df['gene_group'] = filtered_df['gene_id'].astype(str).apply(extract_gene_group)

# Group genes by gene group
summary_rows = []

for group, filtered_group in filtered_df.groupby('gene_group'):
    combined_group = combined_df[combined_df['gene_group'] == group]
    genes = filtered_group['gene_id'].unique()
    if len(genes) != 2:
        continue  # Skip groups that don't have exactly 2 orthologous genes

    gene1, gene2 = genes

    # Starting motif counts
    gene1_starting = combined_group[combined_group['gene_id'] == gene1]['motif_id'].nunique()
    gene2_starting = combined_group[combined_group['gene_id'] == gene2]['motif_id'].nunique()

    # Filtered motif counts
    gene1_filtered = filtered_group[filtered_group['gene_id'] == gene1]['motif_id'].unique()
    gene2_filtered = filtered_group[filtered_group['gene_id'] == gene2]['motif_id'].unique()

    # Shared motifs after filtering
    shared_motifs = set(gene1_filtered) & set(gene2_filtered)

    summary_rows.append({
        'Gene Group': group,
        'Gene 1 ID': gene1,
        'Gene 2 ID': gene2,
        'Starting Motifs (Gene 1)': gene1_starting,
        'Starting Motifs (Gene 2)': gene2_starting,
        'Filtered Motifs (Gene 1)': len(gene1_filtered),
        'Filtered Motifs (Gene 2)': len(gene2_filtered),
        'Shared Motifs After Filter': ', '.join(sorted(shared_motifs)),
        'Shared Motif Count': len(shared_motifs)
    })

# Create final DataFrame
summary_df = pd.DataFrame(summary_rows)

# Save to CSV
summary_df.to_csv('orthologous_gene_motif_summary.csv', index=False)

# Optional: Display table
print(summary_df)
