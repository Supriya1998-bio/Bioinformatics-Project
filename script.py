import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os
import numpy as np

# Ortholog mapping: zebrafish gene (lowercase) -> platyfish gene (lowercase)
ortholog_map = {
    'braf': 'braf',
    'cdkn2a/b': 'cdnk2x',
    'cdkn2aip': 'cdkn2aip',
    'e2f1': 'e2f1',
    'foxo3b': 'foxo3b',
    'mitf': 'mitfa',
    'sox10': 'sox10',
    'stat5a': 'stat5a',
    'tp53bp2': 'tp53bp2b',
    'zeb1a': 'zeb1'
}

# Function to extract gene symbol from sequence_name
def extract_symbol(seq_name: str) -> str:
    parts = seq_name.split("|")
    for part in parts:
        if not part.upper().startswith("ENS"):
            return part.strip().lower()
    return parts[-1].strip().lower()

# Load and filter FIMO files
def load_and_process(fimo_file):
    df = pd.read_csv(fimo_file, sep='\t')
    df = df[df['start'] <= 1000]  # promoter region filter
    df['gene'] = df['sequence_name'].apply(extract_symbol)
    return df

# Load data
zf_file = 'FIMO Zebrafish.tsv'
pf_file = 'FIMO Platyfish.tsv'
zf_df = load_and_process(zf_file)
pf_df = load_and_process(pf_file)

# Filter for significant motifs (q < 0.05)
zf_df_q = zf_df[zf_df['q-value'] < 0.05].copy()
pf_df_q = pf_df[pf_df['q-value'] < 0.05].copy()

# Prepare result container
results = []

# Create folder for venn plots
venn_dir = 'venn_plots'
os.makedirs(venn_dir, exist_ok=True)

# Get gene sets
zf_genes = set(zf_df['gene'])
pf_genes = set(pf_df['gene'])

# Helper function to get TF names from dataframe
def get_tf_names(df):
    tf_set = set()
    for _, row in df.iterrows():
        tf_name = row.get('motif_family') or row.get('motif_alt_id')  # fallback
        if pd.notnull(tf_name):
            tf_set.add(tf_name.strip())
    return tf_set

# --- Enhanced Venn Diagram plotting ---
for zf_gene in sorted(zf_genes):
    if zf_gene not in ortholog_map:
        continue

    pf_gene = ortholog_map[zf_gene]
    if pf_gene not in pf_genes:
        continue

    zf_filtered = zf_df_q[zf_df_q['gene'] == zf_gene]
    pf_filtered = pf_df_q[pf_df_q['gene'] == pf_gene]

    zf_tf_names = get_tf_names(zf_filtered)
    pf_tf_names = get_tf_names(pf_filtered)

    shared_tfs = zf_tf_names.intersection(pf_tf_names)
    shared_count = len(shared_tfs)
    shared_names_str = ", ".join(sorted(shared_tfs))

    results.append({
        'gene': zf_gene,
        'zf_q_filtered_TFs': len(zf_tf_names),
        'pf_q_filtered_TFs': len(pf_tf_names),
        'shared_TF_count': shared_count,
        'shared_TF_names': shared_names_str
    })

    plt.figure(figsize=(6,6))
    venn = venn2([zf_tf_names, pf_tf_names], set_labels=('Zebrafish', 'Platyfish'))

    # Add black circle borders to the venn diagram
    for subset_id in ['10', '01', '11']:
        patch = venn.get_patch_by_id(subset_id)
        if patch is not None:
            patch.set_edgecolor('black')
            patch.set_linewidth(1.5)
            patch.set_alpha(0.7)

    # Style labels
    for label_id in ['10', '01', '11']:
        label = venn.get_label_by_id(label_id)
        if label is not None:
            label.set_fontsize(10)
            label.set_fontweight
            label.set_color('black')

    if shared_tfs:
        shared_text = '\n'.join(sorted(shared_tfs))
        if venn.get_label_by_id('11'):
            venn.get_label_by_id('11').set_text(shared_text)
        else:
            plt.text(0.5, 0.5, shared_text, ha='center', va='center', fontsize=8)

    plt.title(f'Shared TFs for {zf_gene}', fontsize=14, fontweight='bold')
    plt.tight_layout()

    safe_gene_name = zf_gene.replace('/', '_')
    plt.savefig(f"{venn_dir}/{safe_gene_name}_TF_venn.png", dpi=150)
    plt.close()

# Save summary table
summary_df = pd.DataFrame(results)
summary_df.to_csv('melanoma_TF_motif_comparison_summary.tsv', sep='\t', index=False)

print("[+] Saved summary table to melanoma_TF_motif_comparison_summary.tsv")
print(f"[+] Venn diagrams saved to '{venn_dir}' folder.\n")
print("Summary table:")
print(summary_df)


# --- Distance from TSS Analysis ---

# Function to calculate distance from TSS using midpoint
def calc_distance_from_tss(row):
    midpoint = (row['start'] + row['stop']) / 2
    distance = 1000 - midpoint  # since promoter region is 1000bp upstream
    return distance

# Calculate distances
zf_df_q['distance_from_tss'] = zf_df_q.apply(calc_distance_from_tss, axis=1)
pf_df_q['distance_from_tss'] = pf_df_q.apply(calc_distance_from_tss, axis=1)

# Create directory for plots
os.makedirs('distance_plots', exist_ok=True)

# Prepare results container for distance summary
distance_summary_rows = []

# --- Enhanced Bar plots & Dot plots ---
for zf_gene in sorted(zf_df['gene'].unique()):
    if zf_gene not in ortholog_map:
        continue
    pf_gene = ortholog_map[zf_gene]
    if pf_gene not in pf_df['gene'].unique():
        continue

    zf_gene_df = zf_df_q[zf_df_q['gene'] == zf_gene]
    pf_gene_df = pf_df_q[pf_df_q['gene'] == pf_gene]

    # Count motifs by distance groups
    def count_distance_groups(df):
        less_100 = (df['distance_from_tss'] < 100).sum()
        more_100 = (df['distance_from_tss'] >= 100).sum()
        return less_100, more_100

    zf_less_100, zf_more_100 = count_distance_groups(zf_gene_df)
    pf_less_100, pf_more_100 = count_distance_groups(pf_gene_df)

    # Save summary rows
    distance_summary_rows.append({
        'gene': zf_gene,
        'species': 'zebrafish',
        '<100bp': zf_less_100,
        '>=100bp': zf_more_100
    })
    distance_summary_rows.append({
        'gene': pf_gene,
        'species': 'platyfish',
        '<100bp': pf_less_100,
        '>=100bp': pf_more_100
    })

    # --- Bar plot ---
    labels = ['<100bp', '≥100bp']
    counts_zf = [zf_less_100, zf_more_100]
    counts_pf = [pf_less_100, pf_more_100]

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7,5))
    rects1 = ax.bar(x - width/2, counts_zf, width, label='Zebrafish', color='#1f77b4', edgecolor='black')
    rects2 = ax.bar(x + width/2, counts_pf, width, label='Platyfish', color='#ff7f0e', edgecolor='black')

    ax.set_ylabel('Motif Count', fontsize=12, fontweight='bold')
    ax.set_title(f'Motif Counts by Distance from TSS for {zf_gene}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Add value labels on bars
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f'{height}',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontsize=10, fontweight='bold')

    autolabel(rects1)
    autolabel(rects2)

    plt.tight_layout()
    safe_gene_name = zf_gene.replace('/', '_')
    plt.savefig(f'distance_plots/{safe_gene_name}_motif_count_bar_chart.png', dpi=150)
    plt.close()

    # --- Dot plot ---
    fig, axs = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    axs[0].scatter(zf_gene_df['distance_from_tss'], zf_gene_df['q-value'],
                   s=50, c='royalblue', edgecolors='black', alpha=0.75)
    axs[0].set_title(f'Zebrafish: {zf_gene}', fontsize=13, fontweight='bold')
    axs[0].set_xlabel('Distance from TSS (bp)', fontsize=12)
    axs[0].set_ylabel('q-value', fontsize=12)
    axs[0].invert_xaxis()
    axs[0].grid(True, linestyle='--', alpha=0.6)

    axs[1].scatter(pf_gene_df['distance_from_tss'], pf_gene_df['q-value'],
                   s=50, c='darkorange', edgecolors='black', alpha=0.75)
    axs[1].set_title(f'Platyfish: {pf_gene}', fontsize=13, fontweight='bold')
    axs[1].set_xlabel('Distance from TSS (bp)', fontsize=12)
    axs[1].invert_xaxis()
    axs[1].grid(True, linestyle='--', alpha=0.6)

    plt.suptitle(f'Q-value vs Distance from TSS for {zf_gene} / {pf_gene}', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f'distance_plots/{safe_gene_name}_distance_dot_plot.png', dpi=150)
    plt.close()

# Save distance summary table
distance_summary_df = pd.DataFrame(distance_summary_rows)
distance_summary_df.to_csv('distance_summary_appendix.tsv', sep='\t', index=False)

print("[+] Distance plots saved in 'distance_plots' folder.")
print("[+] Distance summary appendix saved to 'distance_summary_appendix.tsv'.")


# --- Protein-DNA Interaction Preparation (Part E) ---

# Collect candidate TF-gene pairs where motifs are <100bp from TSS
interaction_candidates = []

for zf_gene in sorted(zf_df_q['gene'].unique()):
    if zf_gene not in ortholog_map:
        continue
    pf_gene = ortholog_map[zf_gene]
    if pf_gene not in pf_df_q['gene'].unique():
        continue

    # Subset for each gene
    zf_gene_df = zf_df_q[zf_df_q['gene'] == zf_gene]
    pf_gene_df = pf_df_q[pf_df_q['gene'] == pf_gene]

    # Filter <100bp motifs
    zf_close = zf_gene_df[zf_gene_df['distance_from_tss'] < 100]
    pf_close = pf_gene_df[pf_gene_df['distance_from_tss'] < 100]

    # Add to candidates
    for _, row in pd.concat([zf_close, pf_close]).iterrows():
        tf_name = row['motif_alt_id'] if pd.notnull(row['motif_alt_id']) else row['motif_id']
        interaction_candidates.append({
            'gene': row['gene'],
            'species': 'zebrafish' if row['gene'] in zf_genes else 'platyfish',
            'TF': tf_name,
            'motif_id': row['motif_id'],
            'matched_sequence': row['matched_sequence'],
            'distance_from_tss': row['distance_from_tss'],
            'q_value': row['q-value']
        })

# Save candidate list for docking prep
interaction_df = pd.DataFrame(interaction_candidates)
interaction_df.to_csv("protein_dna_interaction_candidates.tsv", sep="\t", index=False)

print("[+] Protein–DNA interaction candidates saved to 'protein_dna_interaction_candidates.tsv'.")
print("[i] Use these TF–gene pairs (<100bp from TSS) for docking with HADDOCK or AutoDock.")

# Keep only the best motif per (gene, TF, species)
best_hits = (
    interaction_df
    .sort_values("q_value")  # sort so lowest q_value comes first
    .groupby(["gene", "TF", "species"], as_index=False)
    .first()  # take the first (best) hit
)

# Save a cleaner file
best_hits.to_csv("protein_dna_best_hits.tsv", sep="\t", index=False)

print("[+] Filtered best motif per TF–gene–species saved to 'protein_dna_best_hits.tsv'")
# Load best_hits table
best_hits = pd.read_csv("protein_dna_best_hits.tsv", sep="\t")

# Function to select top N TFs per gene by lowest q_value
def select_top_TFs(df, top_n=2):
    top_candidates = (
        df.sort_values(['gene', 'q_value'])  # sort by gene then lowest q-value
          .groupby('gene', as_index=False)
          .head(top_n)  # take top N rows per gene
    )
    return top_candidates

# Select top 2 TFs per gene
top_hits = select_top_TFs(best_hits, top_n=2)

# Save table
top_hits.to_csv("protein_dna_top_candidates.tsv", sep="\t", index=False)

print("[+] Top 1–2 TF candidates per gene saved to 'protein_dna_top_candidates.tsv'")
print(top_hits)