## Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: RNA-Seq Pipeline


 RNA-seq workflow **from raw SRA download to trimming with Fastp** on the Negishi HPC cluster.

---
Md Imran Hasan (hasan128) | Yini Li Lab
```markdown
# 🧬 RNA-seq Processing Pipeline — GSE124439 (Hammell 2019)

This repository documents a full end-to-end workflow for processing the GSE124439 RNA-seq dataset  
('PRJNA512012') using the Negishi HPC cluster at Purdue.  
The example focuses on the Frontal Cortex, and Motor Cortex (Medial) sub-regions.

---

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
├        
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
## 📊 Cohort & Subregions (GSE124439)

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
Check Type	Description
Expected SRR	Number of accessions listed in case_srr.txt
Found pairs	FASTQ pairs actually present in the folder
Missing	SRR IDs not yet downloaded or converted
Extra	FASTQs not listed in the metadata file
Unpaired	Single-end files missing _1 or _2 partner
Zero-byte FASTQs	Detects incomplete or failed downloads
---
💡 Run this check after every fasterq-dump batch to confirm integrity before trimming.


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





## I want to auto-detect adapters and trim your FASTQ files without using Slurm, the simplest one-liner approach is to use fastp — it automatically detects adapter sequences, trims low-quality bases, and generates QC reports.

Here’s the one-liner you can run directly inside your fastq/ folder:

mkdir -p ../trim_fastp ../qc_fastp && for f in *_1.fastq; do r=${f%_1.fastq}; fastp -i ${r}_1.fastq -I ${r}_2.fastq -o ../trim_fastp/${r}_1.trimmed.fastq -O ../trim_fastp/${r}_2.trimmed.fastq -h ../qc_fastp/${r}_fastp.html -j ../qc_fastp/${r}_fastp.json --detect_adapter_for_pe -w 8; done



## Here’s the one-liner version that skips already processed files (if both trimmed outputs exist):
```bash
mkdir -p ../trim_fastp ../qc_fastp && for f in *_1.fastq; do r=${f%_1.fastq}; if [[ -f ../trim_fastp/${r}_1.trimmed.fastq && -f ../trim_fastp/${r}_2.trimmed.fastq ]]; then echo "Skipping $r (already trimmed)"; else fastp -i ${r}_1.fastq -I ${r}_2.fastq -o ../trim_fastp/${r}_1.trimmed.fastq -O ../trim_fastp/${r}_2.trimmed.fastq -h ../qc_fastp/${r}_fastp.html -j ../qc_fastp/${r}_fastp.json --detect_adapter_for_pe -w 8; fi; done
```

✅ This will:

Create output folders (if missing)

Skip samples already trimmed

Run fastp only for unprocessed ones



The command that:
✅ auto-detects adapters
✅ trims reads using fastp
✅ saves trimmed FASTQs in trim_fastp/
✅ saves QC reports in qc_fastp/

```bash
fastq/
trim_fastp/
    ├─ SRRXXXXXX_1.trimmed.fastq
    ├─ SRRXXXXXX_2.trimmed.fastq
qc_fastp/
    ├─ SRRXXXXXX_fastp.html
    ├─ SRRXXXXXX_fastp.json
```

## Alignment (Mapping) with HISAT2
Align trimmed paired reads to the hg19 genome using HISAT2 (spliced aligner).
Keep all uniquely mapped, properly paired reads.
Remove duplicates (Picard) and mitochondrial/rRNA reads.

```bash
hisat2-build hg19.fa hg19.fa
```

🧩 The general syntax
```bash
hisat2-build [options] <reference_in> <ht2_base>
```
<reference_in> → the FASTA file containing your genome (e.g. hg19.fa)
<ht2_base> → the base name for the output index files that HISAT2 will create

🧠 So in our case
```bash
hisat2-build hg19.fa hg19.fa
```
means:
Input genome FASTA → hg19.fa
Output index prefix (basename) → also hg19.fa


HISAT2 will then create 8 index files named:
```bash
hg19.fa.1.ht2
hg19.fa.2.ht2
hg19.fa.3.ht2
hg19.fa.4.ht2
hg19.fa.5.ht2
hg19.fa.6.ht2
hg19.fa.7.ht2
hg19.fa.8.ht2

```

Here’s our one-liner HISAT2 alignment command (ready to run inside your trim_fastp/ folder):

```bash
mkdir -p ../hisat2_align && for r1 in *_1.trimmed.fastq; do base=${r1%_1.trimmed.fastq}; echo ">> Aligning $base ..."; hisat2 -p 8 -x /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer/hg19.fa -1 ${base}_1.trimmed.fastq -2 ${base}_2.trimmed.fastq --summary-file ../hisat2_align/${base}_alignment_summary.txt --dta | samtools view -@ 8 -bS -o ../hisat2_align/${base}.bam; done
```
sbatch

```bash
#!/bin/bash
#SBATCH -A yinili               # Your account name
#SBATCH -p cpu                  # Partition to use (adjust based on your cluster)
#SBATCH -N 1                     # Number of nodes
#SBATCH -n 32                    # Number of CPUs (32)
#SBATCH -t 12:00:00              # Time limit (adjust as needed)
#SBATCH -J hisat2_align_case     # Job name
#SBATCH -o hisat2_align_case-%j.out  # Standard output file
#SBATCH -e hisat2_align_case-%j.err  # Standard error file

# Load necessary modules
module load hisat2
module load samtools

# Create directory for output if it doesn't exist
mkdir -p ../hisat2_align

# Start alignment loop
for r1 in *_1.trimmed.fastq; do
  base=${r1%_1.trimmed.fastq}
  echo ">> Aligning $base ..."
  hisat2 -p 32 \
    -x /depot/yinili/data/Li_lab/GSE124439_Hammell2019/Refer/hg19.fa \
    -1 ${base}_1.trimmed.fastq \
    -2 ${base}_2.trimmed.fastq \
    --summary-file ../hisat2_align/${base}_alignment_summary.txt \
    --dta | samtools view -@ 32 -bS -o ../hisat2_align/${base}.bam
done
```

```bash
sbatch hisat2_align.sbatch
```

✅ What it does:
Creates an output folder ../hisat2_align/
Aligns each paired FASTQ (_1/_2) to hg19
Produces:
SRRxxxxxx.bam (aligned reads)
SRRxxxxxx_alignment_summary.txt (alignment stats)



📄 Each summary file contains lines like:

```bash
12345678 reads; of these:
  12000000 (97.2%) were paired; of these:
    11800000 (98.3%) aligned concordantly 0 times
    35000 (0.3%) aligned concordantly exactly 1 time
    30000 (0.25%) aligned concordantly >1 times
98.8% overall alignment rate
```

It reports:
Total reads processed
% paired and mapped concordantly
% discordant alignments
Overall alignment rate


