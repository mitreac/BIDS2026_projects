# Classic Bioinformatics Projects
### Graduate Bioinformatics — Genomic & Transcriptomic Data Analysis

**Course Context:** These projects are designed to build familiarity with major genomic and transcriptomic data repositories (GEO, TCGA) and standard analysis pipelines. Each project should be completable in approximately 3–5 hours. You are encouraged to use University of Michigan-provided AI tools to assist with coding and interpretation (see the [AI Tools section](#ai-assistance-university-of-michigan-resources) at the bottom).

---

## Project Option 1: Meta-Analysis of Respiratory Disease Gene Expression

### Overview

Single-dataset differential expression studies are notoriously sensitive to cohort composition, platform differences, and batch effects. This project asks: *how reproducible are differentially expressed gene (DEG) signatures across independent COVID-19 datasets?* You will perform differential expression analysis on multiple GEO datasets individually and then integrate them into a meta-analysis, critically evaluating where results converge and where they diverge.

---

### Learning Objectives

- Navigate GEO to locate, evaluate, and download RNA-seq count data and metadata
- Perform and critically interpret differential expression analysis
- Understand the effect of batch correction and data integration on results
- Apply gene set and pathway enrichment analysis
- Synthesize multi-dataset results into a biologically coherent interpretation

---

### Required Datasets

Download the following two COVID-19 datasets from GEO. You are **required** to find and add at least one additional COVID-19 GEO dataset of your choice (justify your selection briefly in your report).

| Accession | Description |
|-----------|-------------|
| [GSE202805](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202805) | Required |
| [GSE227116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227116) | Required |
| *(your choice)* | At least 1 additional COVID-19 dataset |

**Data access:** Use the `GEOquery` R package, or download count matrices and metadata directly from the GEO web interface.

---

### Part 1: Data Acquisition and Quality Control

Before any analysis, critically examine each dataset:

- What tissue or cell type was profiled? What sequencing platform was used?
- What are the case and control group definitions? Are they consistent across datasets?
- Perform basic QC: library size distributions, PCA or UMAP colored by dataset, sample, and condition.
- Inspect for and document any obvious batch effects between datasets.

> **Checkpoint question:** Based on QC alone, do you anticipate strong or weak concordance across datasets? Why?

---

### Part 2: Within-Dataset Differential Expression Analysis

Perform differential expression analysis **independently** on each dataset comparing COVID-19 cases to appropriate controls. Use one of the following tools: `DESeq2`, `edgeR`, or `limma/voom`. Apply consistent thresholds across all datasets:

- Adjusted p-value < 0.05
- |log₂ fold change| ≥ 1 (i.e., 2-fold change)

For each dataset, report:

- Total number of DEGs (up- and down-regulated separately)
- A volcano plot
- A heatmap of the top 50 DEGs

---

### Part 3: Meta-Analysis by Dataset Integration

Concatenate the raw count matrices from all datasets into a single combined matrix. Before performing differential expression:

- Apply a batch correction method (e.g., `ComBat-seq` from the `sva` package in R, or `Harmony` in Python) to address inter-dataset technical variation.
- Re-run PCA/UMAP after correction to confirm batch effects are reduced.

Then perform differential expression on the integrated dataset using the same tool and thresholds as Part 2.

**Comparative analysis:**

- Generate a Venn diagram comparing DEG lists across all individual datasets and the meta-analysis result.
- Report the number of DEGs unique to each dataset and the number shared across all datasets.
- Do the meta-analysis results more closely resemble any single dataset? Discuss why.

---

### Part 4: Pathway and Gene Set Enrichment Analysis

Apply gene set enrichment or over-representation analysis to each DEG list (individual datasets and meta-analysis). Acceptable tools include: `fgsea`, `clusterProfiler`, `GSEA`, `SPIA`, or equivalent Python tools (e.g., `GSEApy`).

Use at least one of the following gene set collections: KEGG, Reactome, Hallmark (MSigDB), or GO Biological Process.

- Generate a Venn diagram of dysregulated pathways across datasets and the meta-analysis.
- Report pathways that are consistently identified (across 3+ datasets) and those that are dataset-specific.
- Are the dataset-specific pathways biologically meaningful, or do they appear to reflect noise or cohort-specific confounders?

---

### Part 5: Your Original Research Question

Formulate and answer **one novel question** using the data you have collected. Your question should go beyond replicating the steps above. Examples of the *type* of question expected (do not simply copy these):

- Are the consistently dysregulated genes enriched for any specific biological process not captured by standard pathway databases (e.g., using a custom gene set)?
- Do datasets that share a tissue type show higher DEG concordance than datasets from different tissues?
- Does the direction of fold change for shared DEGs agree across datasets, or are there discordant genes?
- Can the top meta-analysis DEGs distinguish COVID-19 from other respiratory diseases if validated against a publicly available influenza dataset?

Clearly state your question, your analytical approach, your result, and your interpretation.

---

### Deliverables

- A written report (~2–3 pages, not counting figures) covering all five parts
- All figures (volcano plots, heatmaps, Venn diagrams, PCA/UMAP plots)
- A brief methods section describing tools, versions, and parameter choices
- Well-commented code (R script, Python script, or Jupyter/Rmd notebook)

---

---

## Project Option 3: Prostate Cancer Transcriptomic Subtyping

### Overview

Prostate cancer is clinically heterogeneous: patients with similar histological grade can have vastly different outcomes. This project asks: *do molecular subtypes defined from transcriptomic data align with established clinical classifications?* You will download TCGA prostate adenocarcinoma (PRAD) data, perform unsupervised clustering to discover transcriptomic subtypes, evaluate how well those subtypes correspond to clinical labels, and optionally build a supervised classifier.

---

### Learning Objectives

- Access and process TCGA RNA-seq data and associated clinical metadata
- Perform normalization and dimensionality reduction on large expression matrices
- Apply and evaluate unsupervised clustering methods
- Compare data-driven subtypes to clinically defined categories
- Optionally: build and evaluate a supervised machine learning classifier

---

### Part 1: Data Acquisition and Clinical Metadata Curation

Download TCGA-PRAD RNA-seq data (raw counts or FPKM/TPM) using one of the following:

- **R:** `TCGAbiolinks` (preferred) or `TCGA2STAT`
- **Python:** `gdc-client` or `pytcga`

Download the associated clinical metadata. Inspect the metadata to identify relevant clinical grouping variables. You must use **at least two** of the following to define sample groups:

| Variable | Notes |
|----------|-------|
| Gleason score | Standard threshold: ≤6 (low), 7 (intermediate), ≥8 (high) |
| Pathologic T stage | pT2 vs. pT3/T4 |
| PSA level at diagnosis | Common threshold: <10, 10–20, >20 ng/mL |
| Biochemical recurrence status | Available as a binary label |

Document how many samples fall into each group and flag any samples with missing clinical data (decide whether to exclude or impute, and justify your choice).

**If you downloaded raw counts:** normalize to CPM, TPM, or FPKM before proceeding, and briefly explain why normalization is necessary.

---

### Part 2: Dimensionality Reduction and Exploratory Analysis

High-dimensional expression data requires reduction before clustering is interpretable.

- Filter to a biologically informative gene subset (e.g., top 2,000–5,000 most variably expressed genes, or a curated list of known prostate cancer marker genes).
- Apply PCA and either UMAP or t-SNE.
- Visualize samples colored by each of your clinical grouping variables.

> **Checkpoint question:** Do any clinical labels show visible separation in reduced-dimensional space before clustering? What does this suggest about the transcriptomic signal for that variable?

---

### Part 3: Unsupervised Clustering and Subtype Discovery

Apply at least **two** unsupervised clustering approaches from the list below and compare their results:

- k-means clustering
- Hierarchical clustering (with your choice of linkage and distance metric — justify)
- Consensus clustering (e.g., `ConsensusClusterPlus` in R)
- Gaussian mixture models

For k-means or similar methods, use the elbow method, silhouette scores, or gap statistic to select the number of clusters *k*. Do not simply assume that *k* equals the number of clinical groups.

**Evaluation:**

- Compute the Adjusted Rand Index (ARI) or Normalized Mutual Information (NMI) between your cluster assignments and each clinical grouping variable.
- Generate a clustered heatmap of the top differentially expressed genes across clusters.
- Do your data-driven subtypes correspond to a specific clinical variable more than others? Are any clusters clinically ambiguous?

---

### Part 4: Supervised Classification (Required)

Split the dataset into training (70%) and test (30%) sets, stratified by your primary clinical grouping variable. Build a classifier to predict clinical subtype from expression data using at least one of the following:

- **R:** `caret`, `randomForest`, `e1071` (SVM), or `glmnet`
- **Python:** `scikit-learn` (Random Forest, SVM, Logistic Regression, or Elastic Net)

Report on the test set:

- Accuracy, precision, recall, and F1 score (per class and macro-averaged)
- A confusion matrix
- For tree-based or regularized models: feature importance or top contributing genes

> **Reflection question:** Compare the genes most important to your classifier with the genes driving cluster separation in Part 3. Do they overlap? What does this agreement or disagreement tell you?

---

### Part 5: Your Original Research Question

Formulate and answer **one novel question** using the data. Your question must involve an analysis not already performed in Parts 1–4. Examples of the *type* of question expected:

- Are patients whose transcriptomic subtype is discordant with their clinical label (e.g., high Gleason but low-risk cluster) associated with different survival outcomes in the TCGA clinical follow-up data?
- Can you identify a minimal gene signature (e.g., 10–20 genes) that achieves classification performance comparable to using all variable genes?
- Does incorporating somatic mutation burden or copy number variation (also available via TCGAbiolinks) improve subtype prediction?
- How does classifier performance degrade if trained on one clinical variable (e.g., Gleason score) and tested on another (e.g., recurrence status)?

Clearly state your question, analytical approach, result, and interpretation.

---

### Deliverables

- A written report (~2–3 pages, not counting figures) covering all five parts
- All figures (PCA/UMAP, heatmaps, cluster evaluation plots, confusion matrix)
- A brief methods section describing tools, versions, and parameter choices
- Well-commented code (R script, Python script, or Jupyter/Rmd notebook)

---

---

## AI Assistance: University of Michigan Resources

You are **encouraged** to use AI tools to assist with coding, debugging, and interpretation. The following tools are available to U-M students at no additional cost and require no personal subscription:

| Tool | Access | Best for |
|------|--------|----------|
| [Microsoft Copilot](https://copilot.microsoft.com) | Free with U-M email (sign in with uniqname@umich.edu) | Code generation, debugging, explaining error messages |
| [U-M GPT](https://ai.umich.edu/umgpt/) | U-M login required | General analysis planning, interpreting results |

> **Note:** [Perplexity.ai](https://perplexity.ai) and [rtutor.ai](https://rtutor.ai) are also useful free resources, though they are not U-M-specific. You can use them without a subscription for most tasks.

### Optional Agentic AI Component

For students who want to go further, consider using an AI coding assistant in an **agentic** workflow — where the AI iteratively writes, runs, and refines code in response to intermediate results. This is achievable using tools already available to you:

**Option A — GitHub Copilot in VS Code**
U-M students have access to [GitHub Copilot](https://github.com/features/copilot) for free through the [GitHub Student Developer Pack](https://education.github.com/pack). Use it inside VS Code with the Copilot Chat extension to interactively build and debug your analysis pipeline step-by-step.

**Option B — Microsoft Copilot in a Notebook**
Use Copilot to iteratively interpret QC output, suggest next steps, and generate downstream code blocks based on your actual results. Paste figures or summary tables into the chat and ask follow-up analytical questions.

**What makes it "agentic":** Rather than generating all code upfront, you feed real intermediate outputs (PCA plots, error messages, DEG tables) back to the AI and ask it to revise the approach. Document these iterations in your notebook as a mini-log of your AI-assisted workflow.

**Assessment of your AI use:** In your report, include a brief paragraph describing how you used AI tools, what worked well, what required correction, and what you ultimately had to solve yourself.

---

## General Grading Criteria (Both Projects)

| Component | Weight |
|-----------|--------|
| Data acquisition and QC | 15% |
| Core analysis (Parts 2–4) | 40% |
| Original research question (Part 5) | 25% |
| Code quality and reproducibility | 10% |
| Report clarity and interpretation | 10% |
