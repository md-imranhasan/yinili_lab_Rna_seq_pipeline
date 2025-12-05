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
  -o T2T_CHM13v2_gene_expression_counts_multi_fractional_update_hs1_liftoff_genes_12_5_2025.txt \
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
  -f -M --fraction \
  -a /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Repeatmasker/T2T_CHM13v2_hs1_repeatmasker.gtf \
  -F GTF -t exon -g gene_id \
  -o TE_counts_multi_raw_12_5_New.txt \
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
library(DESeq2)

############################################################
## 1. Read featureCounts output
############################################################

# Change this path if needed:
counts <- read.delim(
  "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/T2T_CHM13v2_gene_expression_counts_multi_fractional_update_hs1_liftoff_genes_12_5_2025_filtered_5.txt",
  comment.char = "#",
  check.names = FALSE
)

# Quick sanity check (optional)
# head(counts[, 1:8])

############################################################
## 2. Clean sample column names (remove folder paths)
############################################################

# Columns that are sample BAMs
sample_cols_full <- grep("\\.clean\\.bam$", colnames(counts), value = TRUE)

# Remove directory prefixes like "case/hisat2_t2t_bam/"
clean_names <- basename(sample_cols_full)

# Replace column names
colnames(counts)[match(sample_cols_full, colnames(counts))] <- clean_names

# Re-detect sample columns with cleaned names
sample_cols <- grep("\\.clean\\.bam$", colnames(counts), value = TRUE)

############################################################
## 3. Build count matrix
############################################################

# Some gene IDs can be duplicated → make unique feature IDs
counts$feature_id <- make.unique(as.character(counts$Geneid))

# Extract counts only for sample columns
counts_mat <- as.matrix(counts[, sample_cols])

# Make sure they are integers (DESeq2 requirement)
storage.mode(counts_mat) <- "integer"

# Replace any NA with 0 (should be rare)
counts_mat[is.na(counts_mat)] <- 0L

# Set rownames to unique feature IDs
rownames(counts_mat) <- counts$feature_id

# Optional check:
dim(counts_mat)
counts_mat[1:5, 1:5]

############################################################
## 4. Define case vs control samples
############################################################

# Your sample IDs (exactly as you listed)
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

# Get plain SRR IDs from column names: "SRRxxxxxx.clean.bam" → "SRRxxxxxx"
samples    <- colnames(counts_mat)
sample_ids <- sub("\\.clean\\.bam$", "", samples)

# Initialize condition vector
condition <- rep(NA_character_, length(sample_ids))

# Assign case / control
condition[sample_ids %in% case_ids]    <- "case"
condition[sample_ids %in% control_ids] <- "control"

# Check if any sample was not assigned
if (any(is.na(condition))) {
  warning("Some samples were not assigned to case or control. Check these IDs:\n",
          paste(sample_ids[is.na(condition)], collapse = ", "))
}

# Make condition a factor; control as reference level
condition <- factor(condition, levels = c("control", "case"))

# Build colData
coldata <- data.frame(
  row.names = samples,
  condition = condition
)

# Quick sanity checks
print(table(coldata$condition))
head(coldata)

############################################################
## 5. Create DESeq2 object and run DE analysis
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ condition
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Get results: case vs control
res <- results(dds, contrast = c("condition", "case", "control"))

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Peek at top genes
head(res_ordered)

############################################################
## 6. Save results to CSV
############################################################

write.csv(
  as.data.frame(res_ordered),
  file = "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/T2T_CHM13v2_gene_expression_counts_multi_fractional_update_hs1_liftoff_genes_12_5_2025_DEG.csv"
)

############################################################
## 7. (Optional) VST and PCA plot
############################################################

# Variance stabilizing transform
vsd <- vst(dds, blind = FALSE)

# Basic PCA plot by condition
plotPCA(vsd, intgroup = "condition")


plotMA(res, ylim = c(-5, 5))


# res is your DESeq2 results table
res_df <- as.data.frame(res)

# Basic volcano plot
with(res_df, plot(
  log2FoldChange,
  -log10(pvalue),
  pch = 20,
  col = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "red", "blue"),
  xlab = "log2 Fold Change",
  ylab = "-log10(p-value)",
  main = "Volcano Plot"
))



library(ggplot2)

res_df <- as.data.frame(res)

res_df$threshold <- with(res_df,
                         padj < 0.05 & abs(log2FoldChange) > 1)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = threshold), alpha = 0.6) +
  scale_color_manual(values = c("yellow", "red")) +
  theme_bw() +
  xlab("log2 Fold Change") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot (Simple)") +
  theme(legend.position = "none")


res_df <- subset(res_df, !is.na(pvalue))





dds       # your DESeqDataSet
res       # results(dds, ...)
res_ordered <- res[order(res$padj), ]

# Basic MA plot
plotMA(res, ylim = c(-5, 5), main = "MA plot: case vs control")

res_df <- as.data.frame(res)
res_df <- subset(res_df, !is.na(baseMean) & !is.na(log2FoldChange))

with(res_df, plot(
  log10(baseMean),
  log2FoldChange,
  pch = 20,
  col = "grey50",
  xlab = "log10(mean expression)",
  ylab = "log2 Fold Change (case / control)",
  main = "MA Plot (simple)"
))

# Highlight significant points
with(subset(res_df, padj < 0.05),
     points(log10(baseMean), log2FoldChange, pch = 20, col = "red"))
abline(h = 0, col = "blue", lty = 2)




library(DESeq2)

vsd <- vst(dds, blind = FALSE)

mat <- assay(vsd)
pc  <- prcomp(t(mat))  # PCA on samples

# Percentage variance
var_expl <- (pc$sdev^2) / sum(pc$sdev^2) * 100

# Get condition from colData
cond <- colData(dds)$condition

plot(pc$x[,1], pc$x[,2],
     pch = 19,
     col = ifelse(cond == "case", "red", "blue"),
     xlab = paste0("PC1 (", round(var_expl[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_expl[2], 1), "%)"),
     main = "PCA: samples")
legend("topright", legend = c("case", "control"),
       col = c("red", "blue"), pch = 19)





# Take top 50 DEGs
res_df <- as.data.frame(res_ordered)
res_df <- subset(res_df, !is.na(padj))
topgenes <- rownames(res_df)[1:50]

mat_top <- assay(vsd)[topgenes, ]

# Scale rows (genes) for better contrast
mat_scaled <- t(scale(t(mat_top)))  # z-score per gene

# Simple heatmap
heatmap(mat_scaled,
        Colv = TRUE,
        Rowv = TRUE,
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(50),
        labRow = topgenes,
        main = "Top 50 DEGs (VST, row-scaled)")







# install.packages("pheatmap")   # if not already installed
library(pheatmap)

res_df <- as.data.frame(res_ordered)
res_df <- subset(res_df, !is.na(padj))
topgenes <- rownames(res_df)[1:50]

mat_top <- assay(vsd)[topgenes, ]

# Row-scaled heatmap
pheatmap(mat_top,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 50 DEGs (VST, row-scaled)")





annotation_col <- data.frame(
  condition = colData(vsd)$condition
)
rownames(annotation_col) <- colnames(vsd)

pheatmap(mat_top,
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 50 DEGs (VST, row-scaled)")













#####find out the significant ones################
res_df <- as.data.frame(res)

# Remove NA values (required to avoid errors)
res_df <- subset(res_df, !is.na(padj) & !is.na(log2FoldChange))

# Add gene name column
res_df$gene <- rownames(res_df)


sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
sig



up   <- subset(sig_TEs, log2FoldChange > 1)
down <- subset(sig_TEs, log2FoldChange < -1)

cat("Total significant genes: ", nrow(sig), "\n")
cat("Upregulated genes: ", nrow(up), "\n")
cat("Downregulated genes: ", nrow(down), "\n")

head(sig[order(sig$padj), ])
```





