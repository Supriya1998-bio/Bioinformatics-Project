This repository contains scripts, results, and documentation for the comparative analysis of transcription factor (TF) binding motifs in the 1000 bp upstream promoter regions of melanoma-associated genes in Danio rerio (zebrafish) and Xiphophorus maculatus (platyfish).
The aim is to identify conserved and species-specific regulatory elements that may contribute to differences in melanoma susceptibility.

🛠 Tools & Methods

FIMO (MEME Suite): motif scanning with the JASPAR 2024 vertebrate motif database.

JASPAR 2024: curated transcription factor binding profiles.

MUSCLE + NIMBOSHade: multiple sequence alignment and visualization of promoter regions.

Python (pandas, matplotlib, matplotlib-venn, numpy): motif filtering, statistical analysis, Venn diagrams, bar plots, and scatter plots.

GitHub repository: contains code, summary tables, and figures for reproducibility.

🔬 Workflow Summary

Selected 10 melanoma-associated genes (e.g., mitf, sox10, braf, cdkn2a/b, zeb1a).

Extracted 1000 bp upstream promoter sequences from Ensembl (canonical transcripts).

Aligned orthologous promoter sequences (zebrafish vs platyfish) using MUSCLE and visualized conserved regions in NIMBOSHade.

Scanned promoter regions with FIMO (q-value < 0.05 threshold).

Filtered motifs within the -1000 bp promoter window.

Compared motif sets across species using Venn diagrams.

Analyzed motif proximity to TSS (<100 bp vs ≥100 bp) with bar plots and scatter plots.

Prioritized candidate TF–gene pairs (<100 bp, lowest q-values) for future protein–DNA docking studies.

📊 Key Results

Conservation: Several genes (e.g., zeb1a) showed conserved motifs (BPC1, BPC5, Klf15, PRDM9) across both species.

Species-specific regulation: Some genes (e.g., cdkn2a/b, sox10) displayed minimal overlap, indicating divergent regulation.

Proximity analysis: Platyfish sox10 had a dense cluster of >100 proximal motifs, suggesting strong promoter activity.

Candidate interactions: High-confidence pairs included sox10–BPC1/BPC6 (platyfish), braf–BPC1/BPC5 (zebrafish), and zeb1a–Klf15 (zebrafish).

👉 This repository provides a reproducible framework for comparative promoter motif analysis and identifies candidate transcriptional regulators that may influence melanoma susceptibility in fish models.
