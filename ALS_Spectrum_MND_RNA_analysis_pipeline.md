Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: RNA-Seq Workflow Infographic

Project Overview

This repository contains a single-page interactive infographic visualizing the complex bioinformatics pipeline for identifying distinct molecular subtypes of Amyotrophic Lateral Sclerosis (ALS) using RNA-sequencing data from postmortem cortex samples, with a specialized focus on Transposable Element (TE) analysis.

The infographic is designed to clearly map the process from raw FASTQ data through quality control, parallel quantitation of genes and TEs, differential expression analysis, and final subtype discovery.

Repository Contents

index.html: The complete, single-file, responsive web application containing all HTML, Tailwind CSS styling, and Chart.js/JavaScript logic.

README.md: This file.

.gitignore: Standard file for ignoring common development artifacts.

Viewing the Infographic

Since this is a single HTML file, there is no build process required.

Option 1: Using GitHub Pages (Recommended for Sharing)

Push this repository to GitHub.

Go to Settings > Pages.

Set the source branch to main (or master) and the folder to (root).

Save. GitHub will deploy the site, and the infographic will be accessible via the provided GitHub Pages URL.

Option 2: Local Viewing

Clone or download this repository.

Open the index.html file directly in any modern web browser (e.g., Chrome, Firefox).

Technical Stack

Front-End: Pure HTML5

Styling: Tailwind CSS (via CDN)

Data Visualization: Chart.js (for the mock pathway analysis bar chart)

Design: Responsive (Mobile-first approach) using the "Brilliant Blues" color palette.

Workflow Modules Visualized

Preprocessing & Alignment: Raw data to clean BAM files.

Parallel Quantitation: Separate paths for Standard Genes (FeatureCounts) and Transposable Elements (TEtranscripts).

Discovery & Interpretation: DESeq2 for differential expression, PCA for subtype discovery, and functional analysis (GO/RBP) for biological insight.
