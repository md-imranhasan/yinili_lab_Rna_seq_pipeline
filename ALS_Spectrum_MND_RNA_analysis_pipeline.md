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



### Duplicate Removal Script for Sorted BAMs
```
module load picard/2.27.5  # Load Picard module

# Loop through all sorted BAM files
for f in ../hisat2_t2t_bam/*.sorted.bam; do
  r=$(basename ${f%.sorted.bam})  # Extract sample name

  echo "Marking duplicates for $r"

  # Mark duplicates and remove them
  picard MarkDuplicates \
    I=$f \
    O=../hisat2_t2t_bam/${r}.dedup.bam \
    M=../hisat2_t2t_bam/${r}.dedup.metrics.txt \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Optional: if you want to keep the duplicates marked but not removed
  # REMOVE_DUPLICATES=false
done
```
#### HISAT2 Duplicate Removal and Indexing
This section explains the parameters used for removing duplicates and creating an index in a BAM file after HISAT2 alignment.

- **I**=`$f`: **Input Sorted BAM File**  
  The input BAM file, which is sorted, that contains the aligned reads.

- **O**=`../hisat2_t2t_bam/${r}.dedup.bam`: **Output BAM File (with duplicates removed)**  
  The output file where the duplicates will be removed and the results saved as a new BAM file.

- **M**=`../hisat2_t2t_bam/${r}.dedup.metrics.txt`: **Duplicate Removal Metrics File**  
  A file that logs the metrics related to the duplicate removal process, including how many duplicates were removed.

- **REMOVE_DUPLICATES=true**: **Duplicate Removal Option**  
  When set to `true`, this option ensures that duplicate reads are removed from the BAM file.

- **CREATE_INDEX=true**: **Create BAM Index**  
  This option creates an index for the deduplicated BAM file, allowing for efficient retrieval of alignments.

- **VALIDATION_STRINGENCY=SILENT**: **Validation Stringency**  
  This option silences warnings related to invalid reads, allowing the process to continue without interruptions.



  #### SBATCH Script for Marking (Not Removing) Duplicates

```
#!/bin/bash
#SBATCH -A yang3099                  # Correct account for Gautschi cluster
#SBATCH -p cpu                        # CPU partition
#SBATCH -N 1                           # Use 1 node
#SBATCH -n 16                          # Number of CPU cores (16 for Picard)
#SBATCH -t 12:00:00                    # Time limit for the job (adjust if necessary)
#SBATCH -J picard_mark_duplicates      # Job name
#SBATCH -o picard_mark_duplicates-%j.out  # Output log
#SBATCH -e picard_mark_duplicates-%j.err  # Error log

# Load necessary modules
module --force purge
module load biocontainers
module load picard

# Change to the directory where your sorted BAM files are located
cd "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/motor_cortex_\(lateral\)/case/hisat2_t2t_bam"

# Create directories for output if not already created
mkdir -p ../qc_markdup

# Loop through all sorted BAM files
for f in ../hisat2_t2t_bam/*.sorted.bam; do
  r=$(basename ${f%.sorted.bam})  # Extract sample name from BAM file

  echo "Marking duplicates for $r"

  # Mark duplicates but don't remove them
  picard MarkDuplicates \
    I=$f \
    O=../hisat2_t2t_bam/${r}.dedup.bam \
    M=../qc_markdup/${r}.dedup.metrics.txt \
    REMOVE_DUPLICATES=false \  # Keep the duplicates marked but not removed
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Optional: If you want to keep the duplicates marked but not removed, we already set REMOVE_DUPLICATES=false
done

echo "All duplicates marked but not removed."

```

Correct one-line command for this GTF
```
awk '!/^#/ && $0 ~ /rRNA/ {OFS="\t"; print $1, $4-1, $5, $9}' GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf > chm13_rRNA.bed
```
🔍 check it worked
```
wc -l chm13_rRNA.bed
head chm13_rRNA.bed
```
```
bam="SRR8375275.dedup.bam"; bed="/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/chm13_rRNA.chr.bed"; total=$(samtools view -@4 -c -F2308 "$bam"); mito=$(samtools view -@4 -c -F2308 "$bam" chrM); rrna=$(samtools view -@4 -c -F2308 -L "$bed" "$bam"); mito_pct=$(awk -v m=$mito -v t=$total 'BEGIN{printf "%.4f",(m/t)*100}'); rrna_pct=$(awk -v r=$rrna -v t=$total 'BEGIN{printf "%.4f",(r/t)*100}'); echo -e "Total:\t$total\nMito:\t$mito\t(${mito_pct}%)\nrRNA:\t$rrna\t(${rrna_pct}%)"
```

```bash
cd /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T

awk 'BEGIN{OFS="\t"}
{
  if ($1=="NC_060925.1")      $1="chr1";
  else if ($1=="NC_060926.1") $1="chr2";
  else if ($1=="NC_060927.1") $1="chr3";
  else if ($1=="NC_060928.1") $1="chr4";
  else if ($1=="NC_060929.1") $1="chr5";
  else if ($1=="NC_060930.1") $1="chr6";
  else if ($1=="NC_060931.1") $1="chr7";
  else if ($1=="NC_060932.1") $1="chr8";
  else if ($1=="NC_060933.1") $1="chr9";
  else if ($1=="NC_060934.1") $1="chr10";
  else if ($1=="NC_060935.1") $1="chr11";
  else if ($1=="NC_060936.1") $1="chr12";
  else if ($1=="NC_060937.1") $1="chr13";
  else if ($1=="NC_060938.1") $1="chr14";
  else if ($1=="NC_060939.1") $1="chr15";
  else if ($1=="NC_060940.1") $1="chr16";
  else if ($1=="NC_060941.1") $1="chr17";
  else if ($1=="NC_060942.1") $1="chr18";
  else if ($1=="NC_060943.1") $1="chr19";
  else if ($1=="NC_060944.1") $1="chr20";
  else if ($1=="NC_060945.1") $1="chr21";
  else if ($1=="NC_060946.1") $1="chr22";
  else if ($1=="NC_060947.1") $1="chrX";
  else if ($1=="NC_060948.1") $1="chrY";
  else if ($1=="NC_012920.1") $1="chrM";
  print
}' chm13_rRNA.bed > chm13_rRNA.chr.bed
```










```text

#!/bin/bash
#SBATCH -A yang3099               # Your account name
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J cleanup_correct
#SBATCH -o cleanup_correct-%j.out
#SBATCH -e cleanup_correct-%j.err



module --force purge
module load biocontainers
module load samtools
module load picard

THREADS=16

# Your working directory with sorted BAMs
cd "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/control/hisat2_t2t_bam"

mkdir -p ../qc_markdup

# ✅ FIX THIS TO REAL PATH
RRNA_BED="/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/chm13_rRNA.chr.bed"
MT_CHR="chrM"   # or "MT" if that’s what idxstats shows

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
    # 3) REMOVE rRNA READS USING BED
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
    # 7) PRINT QC + APPEND SUMMARY
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




# Update 11/24/2025 
## Quantify & Annotation repeat RNAs with RepeatMasker (T2T)
#### Prepare a RepeatMasker annotation for counting (SAF for featureCounts)
RepeatMasker is a tool that identifies repetitive DNA elements in a genome. The RepeatMasker annotation is a list of all repetitive elements that have been annotated and categorized in the genome.

Files we have:

Gene GTF (GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf): This file contains the gene-level annotation (genes, transcripts, exons, etc.) for the T2T CHM13 v2.0 genome.

RepeatMasker GTF (T2T_CHM13v2_hs1_repeatmasker.gtf): This file contains the RepeatMasker annotations (TEs like LINEs, LTRs, Simple repeats, etc.).

### 1️⃣ Combine Gene and Repeat Annotations (GTF)


```bash
cd /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/

# Combine both gene and repeat GTFs into one file
cat GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf T2T_CHM13v2_hs1_repeatmasker.gtf > T2T_CHM13v2_combined.gtf
```

Important step: Sort the combined GTF
```bash
sort -k1,1 -k4,4n T2T_CHM13v2_combined.gtf > T2T_CHM13v2_combined_sorted.gtf
```



RepeatMasker GTF or SAF annotation contains:
```text
Element names (e.g., L1MC3#LINE/L1, AluY#SINE/Alu)

Chromosome locations (e.g., chr1, chr22)

Start and end positions of the repetitive elements

Strand information (whether the element is on the + or - strand)
```

```bash
cd /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/Repeatmasker

awk '
BEGIN {
    FS = "\t";
    OFS = "\t";
    print "GeneID","Chr","Start","End","Strand";
}
# Skip comments
/^#/ { next }

{
    chr=$1; start=$4; end=$5; strand=$7;

    # attribute column
    attr=$9;

    # Extract gene_id from the attribute field
    match(attr, /gene_id "([^"]+)"/, arr);
    gene_id = arr[1];

    if (gene_id == "" || gene_id == NULL) {
        gene_id = "TE_"NR;
    }

    # Print SAF
    print gene_id, chr, start, end, strand;
}
' T2T_CHM13v2_hs1_repeatmasker.gtf > T2T_CHM13v2_repeatmasker.saf

```


#### (Optional but recommended) Create FAMILY-level SAF

```bash
awk '
BEGIN {
    FS="\t"; OFS="\t";
    print "GeneID","Chr","Start","End","Strand";
}
!/^#/ {
    chr=$1; start=$4; end=$5; strand=$7;

    # Extract gene_id
    match($9, /gene_id "([^"]+)"/, arr);
    full=arr[1];     # e.g., L1MC3#LINE/L1

    split(full, a, "#");
    fam=a[2];        # e.g., LINE/L1

    if (fam == "" || fam == NULL) fam = "UNKNOWN";

    print fam, chr, start, end, strand;
}
' T2T_CHM13v2_hs1_repeatmasker.gtf > T2T_CHM13v2_repeatmasker_family.saf

```




## featureCounts on RepeatsRNA
Per Family TE count (recommended)

cd /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/case/hisat2_t2t_bam
```bash
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/Repeatmasker/T2T_CHM13v2_repeatmasker.saf \
  -F SAF \
  -o TE_element_counts_with_BC.txt \
  *.clean.bam

```

### For Unique Reads (Exclude Multi-mapped Reads) - Using .SAF
```bash
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/Repeatmasker/T2T_CHM13v2_repeatmasker.saf \
  -F SAF \
  -o TE_unique_counts.txt \
  *.clean.bam
```

### For Fractional Multi-mapping Reads

```bash
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer_T2T/Repeatmasker/T2T_CHM13v2_repeatmasker.saf \
  -F SAF \
  -M --fraction \
  -o TE_fractional_counts.txt \
  *.clean.bam
```

#### For Unique Reads (Exclude Multi-mapped Reads) - Using .gtf
```text
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -o TE_fractional_counts_gft_unique.txt \
  *.clean.bam
```

### For Fractional Multi-mapping Reads
```text
featureCounts -T 16 -p -B -C -s 2 \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -M --fraction 
  -o TE_fractional_counts_gft_unique.txt \
  *.clean.bam
```


Explanation:
```text
-T 16 – 16 threads

-p – paired-end

-s 2 – reverse-stranded, matches your --rna-strandness RF

-F SAF – annotation is in SAF format

*.clean.bam – only use your final cleaned BAMs

If we want - we can use : -M -O --fraction – count reads that map to multiple features and fractionally assign them
```
Output:

TE_family_counts.txt → matrix of TE family counts (rows = family, cols = samples)

