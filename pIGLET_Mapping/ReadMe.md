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
               ├──► Align to GRCz11 or GRCz12tu (full genome BAM)
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

- `dorado` — long read aligner/basecaller
- `samtools` — BAM manipulation
- `seqtk` — FASTQ subsetting

---
 
## 1. Basecalling (The Entry Point)
The pipeline begins with **Dorado**, Oxford Nanopore's high-performance basecaller. This step transforms raw signal data into readable DNA sequences.

*   **Input:** Raw signal data (POD5) directory.
*   **Process:** Dorado utilizes neural networks to interpret the raw electrical current "squiggles" (ionic current changes) and assign the corresponding A, C, G, or T bases.
*   **Output:** Unaligned FASTQ files, which are subsequently merged into a single master file (`full_fastq`) for comprehensive downstream analysis.

### 2. Dual Reference Alignment
The raw Nanopore FASTQ reads are aligned to two targets using `dorado aligner`:
*   **Host Genome:** Standard alignment to the Zebrafish reference (`GRCz11 or GRCz12tu`).
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
  <img width="3547" height="1806" alt="Zfish mapping visual" src="https://github.com/user-attachments/assets/25982d80-ee48-43a6-afe2-d0bab0530872" />
</p>

