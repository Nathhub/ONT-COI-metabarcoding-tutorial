# 🌊 COI Metabarcoding Data Processing – OceanX Around Africa Expedition

Welcome! 👋  
This repository contains the workflow for extracting and processing the **COI (cytochrome oxidase I)** metabarcoding dataset from the eDNA samples collected during the **OceanX Around Africa** expedition and sequenced using an Oxford Nanopore PromethION sequencer.  
The workflow demonstrates how to go from **raw Nanopore PromethION reads** to a final **OTU (Operational Taxonomic Unit)** and **taxonomy table** using an open-source bioinformatics pipeline.

---

## 🧬 About the Data

A total of **18 eDNA samples** were collected by filtering **2 L of seawater** from **Niskin bottles** mounted on a **CTD rosette**.  
Two depths were sampled at each station:
- 🌞 **Surface**
- 🌊 **Deep Chlorophyll Maximum (DCM)**

DNA was extracted using the **QIAGEN DNeasy PowerWater Kit** and amplified by PCR targeting three molecular markers:
- **16S rRNA** → Bacteria and Archaea  
- **18S rRNA** → Protists and other microeukaryotes  
- **COI (Cytochrome Oxidase I)** → Metazoans (animals)  

Barcoded amplicons were pooled per sample using the **ONT Native Barcoding Kit**, and the final libraries were sequenced across **three PromethION flow cells**.

This repository focuses exclusively on the **COI marker** and its downstream processing.

---

## 📘 What You’ll Find Here

### 📄 [`COI_metabarcoding_tutorial.md`](COI_metabarcoding_tutorial.md)
A step-by-step guide explaining how to:
- Combine Nanopore FASTQ files  
- Trim primers and filter low-quality reads  
- Generate consensus sequences (OTUs)  
- Build and normalize an OTU table  
- Assign taxonomy using BLAST  
- Export data for analysis in R  

Each step includes **commands**, **explanations**, and **teaching notes** to help students understand the reasoning behind each part of the workflow.

### 📂 Example Output Files
If provided, these may include:
- `otu_table.primary.tsv` – raw OTU count table  
- `otu_table.primary.cpm.tsv` – normalized OTU abundances (CPM)  
- `taxonomy_table_named.tsv` – BLAST-based taxonomic assignments  
- `all_samples_consensus_nodup.fasta` – unique consensus sequences  

---

## 💡 Next Steps: Analyze in R

Once you have your OTU and taxonomy tables, continue the analysis in **R** using the  
👉 [**phyloseq** package](https://joey711.github.io/phyloseq/)

With phyloseq, you can:
- Merge OTU, taxonomy, and metadata tables  
- Compute alpha and beta diversity  
- Visualize community composition (barplots, ordinations, etc.)  

> ⚙️ **Note:**  
> The bioinformatics processing steps were performed on a **High-Performance Computing (HPC)** system.  
> You can reproduce them on a **Linux terminal** or **WSL**, but be sure to **adjust the number of threads/cores** (`--cores`, `-j`, or `-t` options) according to your computer’s capacity.

---

## 🔗 Tool References

- **[Cutadapt](https://github.com/marcelm/cutadapt)** – Adapter and primer trimming  
- **[Chopper](https://github.com/wdecoster/chopper)** – Quality and length filtering  
- **[Amplicon Sorter](https://github.com/avierstr/amplicon_sorter)** – Amplicon clustering and consensus generation  
- **[Minimap2](https://github.com/lh3/minimap2)** – Read alignment  
- **[Samtools](https://github.com/samtools/samtools)** – Alignment file processing  
- **[MZG COI Database](https://github.com/mbgmbg/MZGdb)** – Reference database for COI taxonomic assignment  

---

## 👩‍🔬 Author

**Nathan [Your Last Name]**  
National Oceanography Centre  
*Created for student training on eDNA metabarcoding workflows*  

🗓️ **Last updated:** November 2025
