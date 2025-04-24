# MIC
This repository contains code and data for identifying and characterizing metastasis-initiating cells (MICs) in ovarian cancer using scRNA-seq.

# 🧬 Regulation of Initiation and Establishment of Ovarian Cancer Metastasis by Metastasis Initiating Cells

This project explores how metastasis-initiating cells (MICs) drive ovarian cancer (OC) metastasis, with a focus on single-cell transcriptomics and cell–cell interactions.

## 🔍 Highlights

- **Model System**: Combined in vitro 3D omental culture and in vivo xenografts to mimic sequential metastasis stages.
- **Data**: Single-cell RNA-seq (10X Genomics) and public datasets (GSE165897, GSE222557, TCGA OV, AOCS OV).
- **MIC Identification**: CC1 subpopulation identified as MICs with high stemness, hybrid EMT status, and poor clinical prognosis.
- **Trajectory Analysis**: Used Monocle3 to map differentiation from MICs to metastatic subpopulations.
- **Transcription Factors**: SCENIC identified trajectory-specific regulators (e.g., SOX4, MYC, KDM5A, RELB).
- **MIC–TME Crosstalk**: CellChat revealed strong interactions between MICs and stromal cells (fibroblasts, mesothelial cells).
- **L-R Interactions**: Key autocrine and paracrine signaling mechanisms identified (e.g., F2RL1, DSG3, ITGAV/ITGB3).
- **Clinical Relevance**: MIC signatures correlate with worse survival in TCGA and AOCS datasets.

## 🧠 Techniques

- scRNA-seq processing: Seurat, MAGIC, inferCNV
- Stemness inference: CytoTRACE, entropy, AUCell, scFEA
- Pseudotime: Monocle3
- Regulatory analysis: SCENIC
- Cell–cell communication: CellChat
- Survival analysis: GSVA, Kaplan-Meier, Cox regression

## 📁 Dataset Accessions

- GSE165897
- GSE222557
- TCGA OV (via Xena)
- AOCS OV (via Xena)

## 📌 Key Insight

MICs orchestrate ovarian cancer metastasis via dynamic transcriptional regulation and niche interactions with stromal cells. Targeting MIC–TME signaling and differentiation pathways offers new therapeutic opportunities.
