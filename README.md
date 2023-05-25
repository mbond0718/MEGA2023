# MEGA2023

Analysis for "Chromatin loop dynamics during cellular differentiation are associated with changes to both anchor and internal regulatory features"
Preprint: https://www.biorxiv.org/content/10.1101/2022.10.31.514600v1

This repository contains the scripts used to analyze the data for this paper in R. The scripts are broken down into: 
* processing: scripts used to process raw data into more processed format (i.e. differential loop calling, differential gene expression, etc)
* analysis: scripts used to process data and produce plots

Output plots are saved to the `plots` directory. 

Raw and processed data for Hi-C (GSE214123GSE213909), RNA-seq (GSE213386), ATAC-seq (GSE213295), and CUT&RUN (GSE213908) are publicly available on GEO and SRA under SuperSeries 213909. 
