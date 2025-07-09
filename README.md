# Bioinformatics-Project
This repository contains scripts, results, and documentation for the comparative analysis of transcription factor (TF) binding motifs in the 1000 bp upstream promoter regions of melanocyte-related and SRY-box–related genes in Danio rerio (zebrafish) and Xiphophorus maculatus (platyfish). The goal is to identify conserved regulatory elements involved in pigmentation and development.

🛠 Tools & Methods
FIMO (MEME Suite): motif scanning using JASPAR 2024 vertebrate database.

JASPAR: curated TF motif database.

MUSCLE + NIMBOSHade: sequence alignment and visualisation.

Python & Biopython: for sequence manipulation and automation.

R/ggplot2 (optional): for plotting filtered motif counts.

Workflow Summary
Extracted 1000 bp upstream promoter regions using NCBI gene coordinates.

Aligned orthologous promoter sequences to assess conservation.

Scanned promoter sequences with FIMO using JASPAR 2024 motifs.

Filtered motifs by q-value (<0.05) and genomic location (within -1000 bp).

Compared motif content across species and gene groups.

Preliminary Results
Alignments show conservation near TSS.

Filtered motifs vary by species and gene group.

Shared motifs identified between orthologous genes, suggesting conserved regulation.
