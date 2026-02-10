# Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: RNA-Seq Pipeline

RNA-seq workflow **from raw SRA download** on the HPC cluster (Purdue University).

## **Md Imran Hasan (hasan128) | Yini Li Lab**

## 🧬 RNA-seq Processing Pipeline — GSE124439 (Hammell 2019)

This repository documents a full end-to-end workflow for processing the GSE124439 RNA-seq dataset  
('PRJNA512012') using the HPC cluster at Purdue.  
This work focuses on the Frontal Cortex, and Motor Cortex (Medial) sub-regions.

This work flowing this paper: 
1. Postmortem cortex samples identify distinct molecular subtypes of ALS: Retrotransposon activation, oxidative stress, and activated glia. https://www.sciencedirect.com/science/article/pii/S221112471931263X?via%3Dihub


## 📁 Directory Layout
## Repository layout
```text

GSE124439_RNAseq/
           # counts per subregion/case-control
├─ metadata/
│  ├─ SraRunTable.csv                # master table (if you keep one)
│  ├─ SraRunTable_Frontal_Cortex.csv
│  ├─ SraRunTable_motor_cortex_lateral.csv
│  ├─ SraRunTable_motor_cortex_medial.csv
│  ├─ Frontal_cortex_case.txt
│  ├─ Frontal_cortex_control.txt
│  ├─ motor_cortex_lateral_case.txt
│  ├─ motor_cortex_lateral_control.txt
│  ├─ motor_cortex_medial_case.txt
│  └─ motor_cortex_medial_control.txt
├─ adapters/
│  └─ TruSeq3-PE.fa                  # adapter file (text; OK to version)
├─ scripts/
│  ├─ download/
│  │  ├─ prefetch_list.sh            # from list -> .sra cache
│  │  └─ fasterq_from_list.sh        # from list -> FASTQ (threads configurable)
│  ├─ qc/
│  │  ├─ summary_check.sh            # your checker; run inside fastq/
│  │  ├─ run_fastqc.sh
│  │  └─ run_multiqc.sh
│  ├─ trim/
│  │  ├─ trimmomatic_params.md       # brief rationale of params
│  │  ├─ run_trimmomatic_case.slurm
│  │  └─ run_trimmomatic_control.slurm
│  └─ util/
│     └─ verify_fastq_vs_metadata.sh # comm-based cross-checks
├─ regions/
│  ├─ Frontal_Cortex/
│  │  ├─ case/
│  │  │  ├─ fastq/                   # (git-ignored)
│  │  │  ├─ trim_trimmomatic/        # (git-ignored)
│  │  │  ├─ qc/                      # (git-ignored)
│  │  │  └─ logs/                    # (git-ignored)
│  │  └─ control/
│  │     ├─ fastq/                   # (git-ignored)
│  │     ├─ trim_trimmomatic/        # (git-ignored)
│  │     ├─ qc/                      # (git-ignored)
│  │     └─ logs/                    # (git-ignored)
│  ├─ motor_cortex_(lateral)/
│  │  ├─ case/ ... (same as above)
│  │  └─ control/ ... 
│  └─ motor_cortex_(medial)/
│     ├─ case/ ... (same as above)
│     └─ control/ ...
                   
```
````
## 📊Cohort & Subregions (GSE124439)

**Total samples: 162**

| Subregion                 | Phenotype                         | Count | Group    |
|--------------------------|-----------------------------------|------:|----------|
| Frontal Cortex           | ALS Spectrum MND                  |   65  | Case     |
| Frontal Cortex           | Non-Neurological Control          |    9  | Control  |
| Motor Cortex (Lateral)   | ALS Spectrum MND                  |   37  | Case     |
| Motor Cortex (Lateral)   | Non-Neurological Control          |    4  | Control  |
| Motor Cortex (Medial)    | ALS Spectrum MND                  |   38  | Case     |
| Motor Cortex (Medial)    | Non-Neurological Control          |    4  | Control  |
| Motor Cortex (unspecified)| ALS Spectrum MND                 |    5  | Case     |

**Totals:** Case = **145**, Control = **17**, Overall = **162**.

> Notes:  
> • “ALS Spectrum MND” is treated as **Case**; “Non-Neurological Control” as **Control**.  
> • The final row (“Motor Cortex – 5 (ALS Spectrum MND)”) is listed as **Motor Cortex (unspecified)** because no lateral/medial label was provided in the summary above.

---

## ⚙️ Environment Setup

Load SRA Toolkit (for downloading) or Trimmomatic / FastQC (for QC).

```bash
module --force purge
module load biocontainers
module load sra-tools/2.11.0-pl5262
# later steps:
module load trimmomatic/0.39
module load fastqc/0.11.9
module load multiqc/1.14
````

---

## 🪄 Step 1 – Configure and Download SRA Data

### 1.1 Prepare directories

```bash
cd /depot/yinili/data/Li_lab/GSE124439_Hammell2019
mkdir -p Frontal_Cortex/{fastq,metadata,qc,trim,logs}
```

### 1.2 Load toolkit & configure

```bash
module load biocontainers
module load sra-tools/2.11.0-pl5262
vdb-config --prefetch-to-cwd   # ensures files download into the current directory
```

### 1.3 Download using accession list

Example: `Frontal_cortex_case.txt` contains one SRR per line.

```bash
while read SRR; do
  echo "Downloading $SRR ..."
  prefetch $SRR
done < Frontal_cortex_case.txt
```

### 1.4 Convert `.sra` → `.fastq`

```bash
mkdir -p fastq
for sra in SRR*/SRR*.sra; do
  echo "Processing $sra ..."
  fasterq-dump --split-files "$sra" -O fastq/
done
```

---

## 🧩 Step 2 – Verify FASTQ Completeness

Inside each `fastq/` folder:

```bash
ls *_1.fastq | sed 's/_1\.fastq//' | sort -u > fastq_srr.txt
awk -F',' 'NR>1{print $1}' ../SraRunTable_Frontal_Cortex.csv | sort -u > metadata_srr.txt
echo "Missing:" && comm -23 metadata_srr.txt fastq_srr.txt
echo "Extra:"   && comm -13 metadata_srr.txt fastq_srr.txt
```
🧾 Step 2.1 – Automatic FASTQ Verification Script
Save this as summary_check.sh inside each fastq/ directory and run it with:
```bash
bash summary_check.sh
```

```bash
# summary_check.sh (run inside the fastq/ folder)
sort -u case_srr.txt > expected.txt
ls -1 *_1.fastq 2>/dev/null | sed 's/_1\.fastq$//' | sort -u > actual_from_fastq.txt

exp=$(wc -l < expected.txt)
found=$(wc -l < actual_from_fastq.txt)
echo "Expected SRR: $exp"
echo "Found pairs : $found"

echo "Missing:"
comm -23 expected.txt actual_from_fastq.txt || true
echo "Extra:"
comm -13 expected.txt actual_from_fastq.txt || true

echo "Unpaired:"
while read SRR; do
  [[ -f ${SRR}_1.fastq && -f ${SRR}_2.fastq ]] || echo "$SRR"
done < expected.txt

echo "Zero-byte FASTQs:"
find . -maxdepth 1 -name "*.fastq" -size 0 -printf "%f\n"
```
🧩 Output Overview

### Check Type Description

- **Expected SRR**  
  Number of accessions listed in `case_srr.txt`.
- **Found pairs**  
  FASTQ pairs actually present in the folder.
- **Missing**  
  SRR IDs not yet downloaded or converted.
- **Extra**  
  FASTQs not listed in the metadata file.
- **Zero-byte FASTQs**  
  Detects incomplete or failed downloads.

> Run this check after every `fasterq-dump` batch to confirm integrity before trimming.


## 🧪 Step 3 – Quality Control (FastQC + MultiQC)

```bash
cd fastq
mkdir -p ../qc
fastqc *.fastq -t 8 -o ../qc
cd ../qc
multiqc .
```

**Interpretation highlights**

| Metric                    | What to check                   | Action                  |
| ------------------------- | ------------------------------- | ----------------------- |
| Per-base quality          | Should stay above Q20–30        | Trim low-quality tails  |
| Adapter content           | Should be near zero             | If high → trim adapters |
| Overrepresented sequences | Detect adapters/primers         | Trim if present         |
| GC content                | Should match expected (~40–50%) | Check for contamination |

---
# 11/13/2025 Update

### I want to auto-detect adapters and trim your FASTQ files without using Slurm, the simplest one-liner approach is to use fastp — it automatically detects adapter sequences, trims low-quality bases, and generates QC reports.

Here’s the one-liner you can run directly inside your fastq/ folder:
```bash
mkdir -p ../trim_fastp ../qc_fastp && for f in *_1.fastq; do r=${f%_1.fastq}; fastp -i ${r}_1.fastq -I ${r}_2.fastq -o ../trim_fastp/${r}_1.trimmed.fastq -O ../trim_fastp/${r}_2.trimmed.fastq -h ../qc_fastp/${r}_fastp.html -j ../qc_fastp/${r}_fastp.json --detect_adapter_for_pe -w 8; done
```


## Here’s the one-liner version that skips already processed files (if both trimmed outputs exist):
```bash
mkdir -p ../trim_fastp ../qc_fastp && for f in *_1.fastq; do r=${f%_1.fastq}; if [[ -f ../trim_fastp/${r}_1.trimmed.fastq && -f ../trim_fastp/${r}_2.trimmed.fastq ]]; then echo "Skipping $r (already trimmed)"; else fastp -i ${r}_1.fastq -I ${r}_2.fastq -o ../trim_fastp/${r}_1.trimmed.fastq -O ../trim_fastp/${r}_2.trimmed.fastq -h ../qc_fastp/${r}_fastp.html -j ../qc_fastp/${r}_fastp.json --detect_adapter_for_pe -w 8; fi; done
```

✅ This will:
```
Create output folders (if missing)
Skip samples already trimmed
Run fastp only for unprocessed ones
```


The command that:
```
✅ auto-detects adapters
✅ trims reads using fastp
✅ saves trimmed FASTQs in trim_fastp/
✅ saves QC reports in qc_fastp/
```
```bash
fastq/
trim_fastp/
    ├─ SRRXXXXXX_1.trimmed.fastq
    ├─ SRRXXXXXX_2.trimmed.fastq
qc_fastp/
    ├─ SRRXXXXXX_fastp.html
    ├─ SRRXXXXXX_fastp.json
```

# Update 11/17/2025
# 🧬 RNA-seq Alignment to T2T-CHM13 Using HISAT2

*(Build Index → Align Reads → Sorted BAM)*

This guide explains how to:

1. Download the **T2T CHM13v2.0 genome**
2. Build the **HISAT2 index**
3. Align **paired-end, strand-specific RNA-seq reads**
4. Produce **sorted BAM** output


---

## 📥 1. Download the T2T Genome

```bash
wget [https://s3.amazonaws.com/nanopore-human-wgs/chm13v2.0.fa.gz](https://s3.amazonaws.com/nanopore-human-wgs/chm13v2.0.fa.gz)
gunzip chm13v2.0.fa.gz
```

This produces:

```
chm13v2.0.fa
```

---

## 🧱 2. Build the HISAT2 Index

```bash
hisat2-build -p 16 chm13v2.0.fa chm13_index
```

This generates eight index files:

```
chm13_index.1.ht2
chm13_index.2.ht2
chm13_index.3.ht2
chm13_index.4.ht2
chm13_index.5.ht2
chm13_index.6.ht2
chm13_index.7.ht2
chm13_index.8.ht2
```

---

## 🎯 3. Align Paired FASTQ Files (Stranded RNA-seq)

Move into trimmed FASTQ directory (e.g., `trim_fastp/`):

```bash
cd trim_fastp/
```

Run the HISAT2 + Samtools pipeline:

```bash
mkdir -p ../hisat2_t2t_bam ../hisat2_t2_tlogs && \
for f in *_1.trimmed.fastq; do
  r=${f%_1.trimmed.fastq}
  echo "Processing: ${r}"
  
  # HISAT2 alignment with paper parameters
  hisat2 -p 16 \
    --rna-strandness RF \
    --dta \
    --summary-file ../hisat2_t2tlogs/${r}.hisat2.summary.txt \
    -x /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/chm13v2.0 \
    -1 ${r}_1.trimmed.fastq \
    -2 ${r}_2.trimmed.fastq \
    2>> ../hisat2_t2tlogs/${r}.hisat2.stderr.log \
  | samtools sort -@ 8 -o ../hisat2_t2t_bam/${r}.sorted.bam
  
  samtools index -@ 8 ../hisat2_t2t_bam/${r}.sorted.bam
  
  # Separate FORWARD strand (flags 99 OR 147)
  samtools view -@ 4 -b -f 99 ../hisat2_t2t_bam/${r}.sorted.bam > ../hisat2_t2t_bam/${r}_f99.bam
  samtools view -@ 4 -b -f 147 ../hisat2_t2t_bam/${r}.sorted.bam > ../hisat2_t2t_bam/${r}_f147.bam
  samtools merge -@ 4 -f ../hisat2_t2t_bam/${r}_forward_strand.bam \
    ../hisat2_t2t_bam/${r}_f99.bam \
    ../hisat2_t2t_bam/${r}_f147.bam
  rm ../hisat2_t2t_bam/${r}_f99.bam ../hisat2_t2t_bam/${r}_f147.bam
  
  # Separate REVERSE strand (flags 83 OR 163)
  samtools view -@ 4 -b -f 83 ../hisat2_t2t_bam/${r}.sorted.bam > ../hisat2_t2t_bam/${r}_f83.bam
  samtools view -@ 4 -b -f 163 ../hisat2_t2t_bam/${r}.sorted.bam > ../hisat2_t2t_bam/${r}_f163.bam
  samtools merge -@ 4 -f ../hisat2_t2t_bam/${r}_reverse_strand.bam \
    ../hisat2_t2t_bam/${r}_f83.bam \
    ../hisat2_t2t_bam/${r}_f163.bam
  rm ../hisat2_t2t_bam/${r}_f83.bam ../hisat2_t2t_bam/${r}_f163.bam
  
  # Index strand-specific BAMs
  samtools index -@ 4 ../hisat2_t2t_bam/${r}_forward_strand.bam
  samtools index -@ 4 ../hisat2_t2t_bam/${r}_reverse_strand.bam
  
  # Quick verification
  echo "  Total: $(samtools view -c ../hisat2_t2t_bam/${r}.sorted.bam)"
  echo "  Forward: $(samtools view -c ../hisat2_t2t_bam/${r}_forward_strand.bam)"
  echo "  Reverse: $(samtools view -c ../hisat2_t2t_bam/${r}_reverse_strand.bam)"
done

```
```text
Output Files

Sorted BAM files: ../hisat2_t2t_bam/ (e.g., SRR8375298.sorted.bam)

Strand-specific BAM files:

Forward strand: SRR8375298_forward_strand.bam

Reverse strand: SRR8375298_reverse_strand.bam

Summary files: ../hisat2_t2tlogs/ (e.g., SRR8375298.hisat2.summary.txt)

```



### BAM Cleanup 
`
# Update 11/21/2025
  #### SBATCH Script 

```
module --force purge
module load biocontainers
module load samtools
module load picard

THREADS=16

# Your working directory with sorted BAMs
cd "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/case/hisat2_t2t_bam"

mkdir -p ../qc_markdup

# ✅ FIX THIS TO REAL PATH
RRNA_BED="/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/chm13_rRNA.chr.bed"
MT_CHR="chrM"   

   # Summary file for all samples
SUMMARY="../qc_markdup/cleanup_summary.tsv"
echo -e "sample\tsorted_reads\tdedup_reads\tduplication_pct\tmt_reads\tmt_pct\trrna_reads\trrna_pct\tclean_reads" > "$SUMMARY"

################################################################################
# MAIN LOOP
################################################################################

for in_bam in *.sorted.bam; do
    sample=${in_bam%.sorted.bam}
    echo "=== Sample: ${sample} ==="

    dedup="${sample}.dedup.bam"
    metrics="../qc_markdup/${sample}.dedup.metrics.txt"
    noMT="${sample}.dedup.noMT.bam"
    clean="${sample}.clean.bam"
    rRNA_only="${sample}.rRNA_only.bam"

    ###########################################################################
    # 1) REMOVE DUPLICATES (MarkDuplicates, NO INDEXING HERE)
    ###########################################################################
    if [[ -f "$dedup" ]]; then
        echo ">> Dedup exists, skipping MarkDuplicates: $dedup"
    else
        echo ">> Running MarkDuplicates on $in_bam ..."
        picard MarkDuplicates \
            -I "${in_bam}" \
            -O "${dedup}" \
            -M "${metrics}" \
            --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES \
            --REMOVE_DUPLICATES true \
            --VALIDATION_STRINGENCY SILENT
    fi

    ###########################################################################
    # 2) REMOVE MITOCHONDRIAL READS (MT / chrM) FROM DEDUP BAM
    ###########################################################################
    if [[ -f "$noMT" ]]; then
        echo ">> noMT exists, skipping MT removal: $noMT"
    else
        echo ">> Removing mitochondrial reads (${MT_CHR} / MT) from $dedup ..."
        samtools view -@ "${THREADS}" -h "${dedup}" \
          | awk -v mt="${MT_CHR}" 'BEGIN{OFS="\t"}
                /^@/ {print; next}
                ($3!=mt && $3!="MT") {print}' \
          | samtools view -@ "${THREADS}" -b -o "${noMT}" -
    fi

    ###########################################################################
    # 3) REMOVE rRNA READS 
    ###########################################################################
    if [[ -f "$clean" && -f "$rRNA_only" ]]; then
        echo ">> Clean and rRNA-only BAMs exist, skipping rRNA removal."
    else
        echo ">> Removing rRNA using BED: ${RRNA_BED}"
        # Reads overlapping RRNA_BED → rRNA_only
        # Reads NOT overlapping (non-rRNA) → clean
        samtools view -@ "${THREADS}" -b -L "${RRNA_BED}" -U "${clean}" "${noMT}" > "${rRNA_only}"
    fi

    ###########################################################################
    # 4) INDEX ONLY FINAL CLEAN BAM
    ###########################################################################
    if [[ ! -f "${clean}.bai" ]]; then
        echo ">> Indexing final clean BAM: ${clean}"
        samtools index -@ "${THREADS}" "${clean}"
    fi

    ###########################################################################
    # 5) QC METRICS — READ COUNTS & PERCENTAGES
    ###########################################################################
    sorted_reads=$(samtools view -@ "${THREADS}" -c "${in_bam}")
    dedup_reads=$(samtools view -@ "${THREADS}" -c "${dedup}")
    noMT_reads=$(samtools view -@ "${THREADS}" -c "${noMT}")
    rrna_reads=$(samtools view -@ "${THREADS}" -c "${rRNA_only}")
    clean_reads=$(samtools view -@ "${THREADS}" -c "${clean}")

    # MT reads removed = dedup - noMT
    mt_reads=$((dedup_reads - noMT_reads))

    # Percentages (protect against division by zero)
    mt_pct=$(awk -v a="$mt_reads" -v b="$dedup_reads" 'BEGIN{if(b==0) print 0; else printf "%.4f", (a/b)*100}')
    rrna_pct=$(awk -v a="$rrna_reads" -v b="$noMT_reads" 'BEGIN{if(b==0) print 0; else printf "%.4f", (a/b)*100}')

    ###########################################################################
    # 6) DUPLICATION % FROM PICARD METRICS (PERCENT_DUPLICATION COLUMN)
    ###########################################################################
    dup_pct=$(awk '
      BEGIN{FS="\t"; col=-1}
      $1=="LIBRARY"{
          for(i=1;i<=NF;i++){
              if($i=="PERCENT_DUPLICATION"){ col=i }
          }
      }
      $1!~/^#|LIBRARY/ && col>0 {
          print $col; exit
      }
    ' "${metrics}")

    ###########################################################################
    # 7) PRINT QC + SUMMARY
    ###########################################################################
    echo "sorted_reads:  ${sorted_reads}"
    echo "dedup_reads:   ${dedup_reads}"
    echo "duplication %: ${dup_pct}"
    echo "mt_reads:      ${mt_reads}  (${mt_pct}%)"
    echo "rRNA_reads:    ${rrna_reads} (${rrna_pct}%)"
    echo "clean_reads:   ${clean_reads}"
    echo

    echo -e "${sample}\t${sorted_reads}\t${dedup_reads}\t${dup_pct}\t${mt_reads}\t${mt_pct}\t${rrna_reads}\t${rrna_pct}\t${clean_reads}" >> "$SUMMARY"

done

echo "Cleanup complete."
echo "Summary table saved to: ${SUMMARY}"
```

### For each sample, the script:

### 📝 Output Summary

For each sample, the pipeline:

- **Appends a tab-separated row** to `cleanup_summary.tsv` containing the following columns:

| Column Name        | Description |
|--------------------|-------------|
| **sample**         | Sample name (derived from `SAMPLE.sorted.bam`) |
| **sorted_reads**   | Number of reads in the original sorted BAM |
| **dedup_reads**    | Number of reads after duplicate removal |
| **duplication_pct** | Duplicate percentage from Picard (`PERCENT_DUPLICATION`) |
| **mt_reads**       | Number of mitochondrial reads removed |
| **mt_pct**         | Percentage of mitochondrial reads (`mt_reads / dedup_reads * 100`) |
| **rRNA_reads**     | Number of rRNA reads removed |
| **rRNA_pct**       | Percentage of rRNA reads  |
| **clean_reads**    | Final number of cleaned reads (no duplicates, no MT, no rRNA) |




### Sample Table:

| Sample      | Sorted Reads | Dedup Reads | Duplication % | MT Reads | MT %   | rRNA Reads | rRNA %  | Clean Reads |
|-------------|--------------|-------------|----------------|----------|--------|------------|---------|-------------|
| SRR8375282  | 78,496,021   | 72,256,630  | 0.089284       | 1,121,960 | 1.5527 | 731,516    | 1.0284  | 70,403,154  |
| SRR8375310  | 100,451,524  | 92,749,156  | 0.085320       | 1,238,321 | 1.3351 | 630,016    | 0.6885  | 90,880,819  |
| SRR8375326  | 98,123,332   | 88,905,654  | 0.101420       | 1,255,705 | 1.4124 | 684,128    | 0.7805  | 86,965,821  |
| SRR8375369  | 61,808,627   | 53,276,718  | 0.160193       | 1,781,408 | 3.3437 | 102,436    | 0.1989  | 51,392,874  |
| SRR8375373  | 82,440,788   | 75,317,170  | 0.099661       | 1,503,821 | 1.9967 | 252,824    | 0.3425  | 73,560,525  |
| SRR8375375  | 80,545,865   | 72,939,064  | 0.106941       | 1,787,054 | 2.4501 | 817,323    | 1.1487  | 70,334,687  |
| SRR8375381  | 109,291,392  | 101,880,994 | 0.075400       | 1,157,733 | 1.1364 | 183,548    | 0.1822  | 100,539,713 |




# Update 12/5/2025

## Quantify & Annotation repeat RNAs with RepeatMasker (T2T)

 ## For genes we are using T2T_CHM13v2_hs1_liftoff_genes.gtf
### For Fractional Multi-mapping Reads
```bash
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/T2T_CHM13v2_hs1_liftoff_genes.gtf \
  -F GTF -t exon -g gene_id \
  -f -M --fraction \
  -o T2T_CHM13v2_gene_expression_counts_multi_fractional_update_hs1_liftoff_genes_12_5_2025.txt \
  case/hisat2_t2t_bam/*.clean.bam \
  control/hisat2_t2t_bam/*.clean.bam
```

### For Unique Reads

```bash
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/T2T_CHM13v2_hs1_liftoff_genes.gtf \
  -F GTF -t exon -g gene_id \
  -f \
  -o T2T_CHM13v2_gene_expression_counts_Unique_update_hs1_liftoff_genes_12_5_2025.txt \
  case/hisat2_t2t_bam/*.clean.bam \
  control/hisat2_t2t_bam/*.clean.bam
```




#### Prepare a RepeatMasker annotation for counting 
RepeatMasker is a tool that identifies repetitive DNA elements in a genome. The RepeatMasker annotation is a list of all repetitive elements that have been annotated and categorized in the genome.
# TE Counting (Fractional Multi-mapping Reads) (corrected 12_5_2025)
### For Fractional Multi-mapping Reads

```text
featureCounts -T 16 -p -B -C -s 2 \
  -f -M --fraction \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -o TE_counts_multi_raw_12_5_New.txt \
  case/hisat2_t2t_bam/*.clean.bam control/hisat2_t2t_bam/*.clean.bam
```

### For Unique Reads

```bash
featureCounts -T 16 -p -B -C -s 2 \
  -f \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -o TE_counts_Unique_raw_12_5_New.txt \
  case/hisat2_t2t_bam/*.clean.bam control/hisat2_t2t_bam/*.clean.bam
```


## Filter the gene count file

```bash
# 1. Load the data
file_path <- "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/T2T_CHM13v2_gene_expression_counts_multi_fractional_update_hs1_liftoff_genes_12_5_2025.txt"

gene_counts <- read.delim(
  file_path,
  comment.char = "#",
  check.names = FALSE
)

# 2. Extract ONLY the count columns (Column 7 to the end)
#    Columns 1-6 are metadata (Geneid, Chr, etc.)
count_cols_index <- 7:ncol(gene_counts)
count_mat <- as.matrix(gene_counts[, count_cols_index])

# 3. Apply the Strict Filter
#    Logic: Count how many samples have >= 5 reads.
#    Keep IF that number equals the TOTAL number of sample columns.
min_reads <- 5
total_samples <- ncol(count_mat)

# Create boolean matrix (TRUE/FALSE)
pass_matrix <- count_mat >= min_reads

# Sum the TRUEs per row
samples_passing <- rowSums(pass_matrix)

# Check if ALL samples passed
keep_genes <- samples_passing == total_samples

# 4. Print Summary
cat("Original Gene Count:", nrow(gene_counts), "\n")
cat("Genes Kept (>=5 reads in EVERY sample):", sum(keep_genes), "\n")
cat("Genes Removed:", sum(!keep_genes), "\n")

# 5. Subset the original dataframe
#    This keeps columns 1-6 intact and filters the rows.
gene_counts_filtered <- gene_counts[keep_genes, ]

# 6. Save to file
write.table(
  gene_counts_filtered,
  file = "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/T2T_CHM13v2_gene_expression_counts_multi_fractional_update_hs1_liftoff_genes_12_5_2025_filtered_5.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```



# Planning to DO NEXT STEPS
#### Gene Expression and TE-Level Counting: Downstream Analysis
### Quality Control (QC) of Count Data


```bash
################################################################################
## RNA-seq Analysis Pipeline: DESeq2
################################################################################

# 0. Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

############################################################
## 1. Read input data
############################################################

# Load the file
counts <- read.delim(
  "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/T2T_CHM13v2_gene_expression_counts_fractional_update_hs1_liftoff_genes_12_5_2025_filtered_5.txt",
  comment.char = "#",
  check.names = FALSE
)

############################################################
## 2. Clean sample column names
############################################################

# Identify columns containing count data (ending in .clean.bam)
sample_cols_full <- grep("\\.clean\\.bam$", colnames(counts), value = TRUE)

# Clean the names (remove directory paths)
clean_names <- basename(sample_cols_full)
colnames(counts)[match(sample_cols_full, colnames(counts))] <- clean_names

# Re-detect sample columns using the new clean names
sample_cols <- grep("\\.clean\\.bam$", colnames(counts), value = TRUE)

############################################################
## 3. Build count matrix
############################################################

# Ensure feature IDs are unique
counts$feature_id <- make.unique(as.character(counts$Geneid))

# Extract count data matrix
counts_mat <- as.matrix(counts[, sample_cols])

# CRITICAL: DESeq2 requires integers. 
# Since input is fractional, we assume integer mode (truncates/rounds implicitly)
storage.mode(counts_mat) <- "integer"

# Replace NAs with 0
counts_mat[is.na(counts_mat)] <- 0L

# Set rownames
rownames(counts_mat) <- counts$feature_id

cat("Initial Matrix Dimensions:", dim(counts_mat), "\n")

############################################################
## 4. Define Case vs Control & Subset Matrix (CRITICAL FIX)
############################################################

# Define your sample lists
case_ids <- c(
  "SRR8375274","SRR8375275","SRR8375276","SRR8375277","SRR8375278",
  "SRR8375279","SRR8375280","SRR8375281","SRR8375283","SRR8375284",
  "SRR8375285","SRR8375286","SRR8375307","SRR8375308","SRR8375309",
  "SRR8375311","SRR8375312","SRR8375324","SRR8375325","SRR8375327",
  "SRR8375328","SRR8375329","SRR8375331","SRR8375333","SRR8375335",
  "SRR8375336","SRR8375337","SRR8375338","SRR8375341","SRR8375345",
  "SRR8375347","SRR8375357","SRR8375362","SRR8375363","SRR8375364",
  "SRR8375365","SRR8375366","SRR8375370","SRR8375371","SRR8375372",
  "SRR8375374","SRR8375376","SRR8375383","SRR8375385","SRR8375386",
  "SRR8375387","SRR8375388","SRR8375389","SRR8375390","SRR8375391",
  "SRR8375392","SRR8375393","SRR8375394","SRR8375395","SRR8375396",
  "SRR8375411","SRR8375414","SRR8375418","SRR8375421","SRR8375424",
  "SRR8375429","SRR8375432","SRR8375437","SRR8375445","SRR8375448"
)

control_ids <- c(
  "SRR8375282","SRR8375310","SRR8375326","SRR8375369","SRR8375373",
  "SRR8375381","SRR8375382","SRR8375384","SRR8375375"
)

# Extract plain sample IDs from columns
current_cols <- colnames(counts_mat)
sample_ids_stripped <- sub("\\.clean\\.bam$", "", current_cols)

# Create metadata frame
coldata <- data.frame(
  row.names = current_cols,
  SampleID  = sample_ids_stripped,
  condition = NA # Initialize as NA
)

# Assign conditions
coldata$condition[coldata$SampleID %in% case_ids]    <- "case"
coldata$condition[coldata$SampleID %in% control_ids] <- "control"

# --- THE FIX STARTS HERE ---
# Remove samples from the matrix that are not in your case/control lists
# (This prevents NAs in the design, which crashes DESeq2)
keep_samples <- !is.na(coldata$condition)

counts_mat <- counts_mat[, keep_samples]
coldata    <- coldata[keep_samples, ]

# Make condition a factor and SET REFERENCE LEVEL
# This ensures "control" is the baseline (Denominator)
coldata$condition <- factor(coldata$condition, levels = c("control", "case"))
coldata$condition <- relevel(coldata$condition, ref = "control")

cat("Final Dimensions for Analysis:", dim(counts_mat), "\n")
print(table(coldata$condition))
# --- THE FIX ENDS HERE ---

############################################################
## 5. Run DESeq2
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ condition
)

# Run the pipeline
dds <- DESeq(dds)

# Get Results
res <- results(dds) # Automatically uses Case vs Control because we set the reference

# Order by significance
res_ordered <- res[order(res$padj), ]

# Export Results
write.csv(as.data.frame(res_ordered), 
          file = "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/T2T_CHM13v2_gene_counts_DESeq2_results_all_filtered10_12_3.csv")

############################################################
## 6. Visualization
############################################################

# --- A. PCA Plot ---
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot: Case vs Control") +
  theme_bw()

# --- B. Volcano Plot ---
res_df <- as.data.frame(res)
# Create a logical column for coloring (padj < 0.05 and FoldChange > 2 aka log2FC > 1)
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")
res_df$Significant[is.na(res_df$Significant)] <- "No"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggtitle("Volcano Plot")

# --- C. Heatmap of Top 50 Genes ---
res_clean <- subset(res_ordered, !is.na(padj))
top_genes <- rownames(res_clean)[1:50]
mat_top   <- assay(vsd)[top_genes, ]

anno_col <- data.frame(Condition = coldata$condition)
rownames(anno_col) <- colnames(mat_top)

pheatmap(mat_top,
         scale = "row",
         annotation_col = anno_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 50 DEGs (VST Z-score)")

############################################################
## 7. Summary Counts
############################################################

sig_res <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
up_genes   <- subset(sig_res, log2FoldChange > 1)
down_genes <- subset(sig_res, log2FoldChange < -1)

cat("----------------------------------\n")
cat("Summary of Significant Genes (padj < 0.05, |LFC| > 1):\n")
cat("Total Significant:", nrow(sig_res), "\n")
cat("Upregulated:      ", nrow(up_genes), "\n")
cat("Downregulated:    ", nrow(down_genes), "\n")
cat("----------------------------------\n")


```

# Update 12/10/2025 (TE corrected code with Filtering )

```bash
library(DESeq2)
library(ggplot2)
library(pheatmap)

te_path <- "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/motor_cortex_(lateral)/TE_counts_unique_raw_12_5_New.txt"

# Skip the first line that starts with "# Program: featureCounts..."
# IMPORTANT: comment.char = "" so that "#" inside Geneid is kept
te <- read.delim(
  te_path,
  skip        = 1,        # skip only the program line
  header      = TRUE,
  check.names = FALSE,
  comment.char = ""       
)

cat("TE table dimensions:", dim(te), "\n")
str(te[1:5, 1:10])









# Columns 1–6 are metadata, 7+ are counts
data_cols <- 7:ncol(te)
sample_cols_full <- colnames(te)[data_cols]

counts_mat <- as.matrix(te[, data_cols])

# Make sure it is numeric
mode(counts_mat) <- "numeric"
counts_mat[is.na(counts_mat)] <- 0

# Row / col names
colnames(counts_mat) <- basename(sample_cols_full)
colnames(counts_mat) <- sub("\\.clean\\.bam$", "", colnames(counts_mat))
rownames(counts_mat) <- make.unique(as.character(te$Geneid))

# Diagnostics
cat("Matrix dimensions (features x samples):", dim(counts_mat), "\n")
cat("Total reads:", sum(counts_mat), "\n")
cat("Non-zero rows:", sum(rowSums(counts_mat) > 0), "\n")
cat("Non-zero columns:", sum(colSums(counts_mat) > 0), "\n")


counts_mat[counts_mat < 5] <- 0


cat("Total reads after <5→0:", sum(counts_mat), "\n")
cat("Non-zero rows after filter:", sum(rowSums(counts_mat) > 0), "\n")




# Reuse original full names for condition
conditions_raw <- ifelse(grepl("^case/", sample_cols_full), "case",
                         ifelse(grepl("^control/", sample_cols_full), "control", NA))
sample_ids <- sub("\\.clean\\.bam$", "", basename(sample_cols_full))

coldata <- data.frame(
  row.names = sample_ids,
  SampleID  = sample_ids,
  condition = conditions_raw,
  stringsAsFactors = FALSE
)

# Align to counts_mat columns
coldata <- coldata[colnames(counts_mat), , drop = FALSE]

# Drop samples with NA condition (if any)
coldata <- coldata[!is.na(coldata$condition), , drop = FALSE]
counts_mat <- counts_mat[, rownames(coldata), drop = FALSE]

cat("Condition table:\n")
print(table(coldata$condition))

coldata$condition <- factor(coldata$condition, levels = c("control", "case"))

# Remove all-zero rows
counts_mat <- counts_mat[rowSums(counts_mat) > 0, , drop = FALSE]
cat("TEs remaining after all-zero row removal:", nrow(counts_mat), "\n")

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ condition
)

cat("Column sums in DESeqDataSet:\n")
print(colSums(counts(dds)))

dds <- DESeq(dds)
res <- results(dds)
res_ordered <- res[order(res$padj), ]

write.csv(as.data.frame(res_ordered),
          file = "TE_DESeq2_motor_lateral_case_vs_control.csv")

```
