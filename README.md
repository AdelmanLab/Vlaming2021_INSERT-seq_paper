# Vlaming2021_INSERT-seq_paper

A collection of scripts that were run to analyse INSERT-seq data and create (tables for) all figures in the Vlaming 2021 manuscript (currently available at https://www.biorxiv.org/content/10.1101/2021.06.01.446655v1)
Scripts are provided for transparency purposes, but path names have not been updated to run for everyone immediately. So, if you want to re-run analyses yourself, you may need to change some path names.

Used in the following order:
1. Mapping scripts (start with WT_uaRNAscreens_sortseq-totalRNA-nascentRNA_mapping).
2. MOODS script to call 5'SS motif matches in inserts.
3. Outputs from mapping used in R analysis 'basis' scripts.
4. Export_tablesforplots_and_GEO-supplements.Rmd, which created all tables used to make plots.
5. Homer script to find enriched motifs in top-10% of TSS-proximal mRNA regions.

Also provided is a short script to generate browser snapshots using the Fluff package.
