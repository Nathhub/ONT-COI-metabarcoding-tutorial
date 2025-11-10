# 🌊 COI Metabarcoding Data Processing – OceanX Around Africa Expedition

Welcome! 👋  
This repository contains the workflow for extracting and processing the **COI (cytochrome oxidase I)** metabarcoding dataset from the eDNA samples collected during the **OceanX Around Africa** expedition and sequenced using an Oxford Nanopore PromethION sequencer.  
The workflow demonstrates how to go from **raw Nanopore PromethION reads** to a final **OTU (Operational Taxonomic Unit)** and **taxonomy table** using a series of open-source bioinformatics tools.

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

---

## 👩‍🔬 Author

**Nathan Hubot**  
National Oceanography Centre  
*Created for student training on eDNA metabarcoding using ONT technology*  

🗓️ **Last updated:** November 2025
