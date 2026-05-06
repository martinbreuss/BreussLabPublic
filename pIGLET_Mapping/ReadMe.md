# Automated mapping of random transgenic insertions in zebrafish for the development of pIGLET lines


## Pipeline Logic: Transgene Integration Site Mapping

This pipeline identifies the exact genomic coordinates of transgene insertions by isolating chimeric reads—those that contain both transgene and host genomic DNA.

## Overview

```
POD5 Files
     │
     └──► Dorado Basecaller (raw signal → DNA sequences)
               │
               Raw FASTQ/BAM
               │
               ├──► Align to GRCz11 (full genome BAM)
               │
               └──► Align to transgene FASTA
                         │
                         └──► mpileup at transgene ends
                                   │
                                   └──► Extract read names at transgene ends
                                             │
                                             └──► Pull those reads from full FASTQ
                                                       │
                                                       └──► Re-align to full genome
                                                                 │
                                                                 └──► Insertion site(s)
```

---

## Requirements

- `dorado` v1.2.0 — long read aligner/basecaller
- `samtools` v1.6 — BAM manipulation
- `seqtk` v1.2 — FASTQ subsetting

---
 
## 1. Basecalling (The Entry Point)
The pipeline begins with **Dorado**, Oxford Nanopore's high-performance basecaller. This step transforms raw signal data into readable DNA sequences.

*   **Input:** Raw signal data (POD5) directory.
*   **Process:** Dorado utilizes neural networks to interpret the raw electrical current "squiggles" (ionic current changes) and assign the corresponding A, C, G, or T bases.
*   **Output:** Unaligned FASTQ files, which are subsequently merged into a single master file (`full_fastq`) for comprehensive downstream analysis.

### 2. Dual Reference Alignment
The raw Nanopore FASTQ reads are aligned to two targets using `dorado aligner`:
*   **Host Genome:** Standard alignment to the Zebrafish reference (`GRCz11`).
*   **Transgene Sequence:** Targeted alignment to the `Transgene` insert sequence.

### 3. Isolating Reads at Transgene boundary with Mpileup
To find the insertion point, the pipeline "baits" reads overlapping the transgene boundaries:
*   Runs `samtools mpileup` on the transgene-aligned BAM.
*   Uses `-Q 0 -q 0` to ensure soft-clipped segments (the parts hanging off into the host genome) are not filtered out.
*   Extracts the **first 25** and **last 25** lines of the mpileup, targeting the start and end of the transgene sequence.

### 4. Read Recovery and Extraction
*   Parses the Read IDs from Transgene-Reference junction points.
*   Uses `seqtk subseq` to retrieve the full-length original sequences for these specific junction-spanning reads from the master FASTQ.

### 5. Integration Localization
*   The isolated junction reads are re-mapped back to the **Host Reference Genome**.
*   The final output is a coordinate-sorted SAM file containing only the reads that bridge the transgene and the genome.

### 6. Final Analysis
The insertion site is identified by viewing the output:
Clusters of these specific reads indicate the precise chromosomal location where the transgene incorporated into the host DNA.

<p align="center">
  <img width="818" height="416" alt="Screenshot 2026-05-06 at 11 03 48 AM" src="https://github.com/user-attachments/assets/dcb2a45d-bd3d-41a4-b3f8-05c7444170e9" />
</p>

