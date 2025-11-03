Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: RNA-Seq Pipeline


 RNA-seq workflow **from raw SRA download to trimming with Trimmomatic** on the Negishi HPC cluster.



---

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
  fasterq-dump "$sra" -O fastq/
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

## ✂️ Step 4 – Trimming (Trimmomatic)

### 4.1 Adapter location

```
/depot/yinili/data/Li_lab/GSE124439_Hammell2019/adapters/TruSeq3-PE.fa
```

### 4.2 Batch trimming script (`run_trimmomatic_case.slurm`)

```bash
#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J trim_case
#SBATCH -o trim_case-%j.out
#SBATCH -e trim_case-%j.err

module --force purge
module load biocontainers
module load trimmomatic/0.39

cd "/depot/yinili/data/Li_lab/GSE124439_Hammell2019/motor_cortex_(medial)/case/fastq"
mkdir -p ../trim_trimmomatic ../logs

for r1 in *_1.fastq; do
  base=${r1%_1.fastq}
  if [[ -f "../trim_trimmomatic/${base}_1.paired.fq.gz" && -f "../trim_trimmomatic/${base}_2.paired.fq.gz" ]]; then
      echo ">> Skipping $base (already trimmed)"
  else
      echo ">> Trimming $base ..."
      trimmomatic PE -threads 8 -phred33 \
        "${base}_1.fastq" "${base}_2.fastq" \
        "../trim_trimmomatic/${base}_1.paired.fq.gz" "../trim_trimmomatic/${base}_1.unpaired.fq.gz" \
        "../trim_trimmomatic/${base}_2.paired.fq.gz" "../trim_trimmomatic/${base}_2.unpaired.fq.gz" \
        ILLUMINACLIP:/depot/yinili/data/Li_lab/GSE124439_Hammell2019/adapters/TruSeq3-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36 \
        2> "../logs/${base}.trimmomatic.log"
  fi
done
```

Submit:

```bash
sbatch run_trimmomatic_case.slurm
```

---

Here’s the **improved one-liner** version:

```bash
for r1 in *_1.fastq; do base=${r1%_1.fastq}; if [[ -f "../trim_trimmomatic/${base}_1.paired.fq.gz" && -f "../trim_trimmomatic/${base}_2.paired.fq.gz" ]]; then echo ">> Skipping $base (already trimmed)"; else echo ">> Trimming $base ..."; trimmomatic PE -threads 8 -phred33 "${base}_1.fastq" "${base}_2.fastq" "../trim_trimmomatic/${base}_1.paired.fq.gz" "../trim_trimmomatic/${base}_1.unpaired.fq.gz" "../trim_trimmomatic/${base}_2.paired.fq.gz" "../trim_trimmomatic/${base}_2.unpaired.fq.gz" ILLUMINACLIP:/depot/yinili/data/Li_lab/GSE124439_Hammell2019/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36 2> "../logs/${base}.trimmomatic.log"; fi; done
```

---

### 🧩 What This Does

* Loops through all `*_1.fastq` files.
* Checks if **both** trimmed paired files exist:

  ```
  ../trim_trimmomatic/${base}_1.paired.fq.gz
  ../trim_trimmomatic/${base}_2.paired.fq.gz
  ```
* ✅ If both exist → **skip** the sample and print:

  ```
  >> Skipping SRR837xxxx (already trimmed)
  ```
* 🚀 Otherwise → runs **Trimmomatic** and logs the output.

---

## 🧭 Job Monitoring on Negishi

```bash
squeue -u $USER
squeue -j <jobid> -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,AllocCPUS
```

* **PD (Priority)** → waiting in queue
* **R (Running)** → active
* **PartitionDown** → cluster maintenance
* **Resources** → waiting for free nodes

---

## 🧠 Trimming Parameters Explained

| Parameter                            | Purpose                                   |
| ------------------------------------ | ----------------------------------------- |
| `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10` | Removes Illumina adapters (≤2 mismatches) |
| `SLIDINGWINDOW:4:20`                 | Cuts when 4-bp window drops below Q20     |
| `LEADING:3 TRAILING:3`               | Trims low-quality ends                    |
| `MINLEN:36`                          | Discards reads shorter than 36 bp         |

---

## ✅ Step 5 – Post-Trim Quality Check (Optional)

```bash
module load fastqc/0.11.9 multiqc/1.14
mkdir -p ../qc_trim
fastqc ../trim_trimmomatic/*paired.fq.gz -t 8 -o ../qc_trim
cd ../qc_trim && multiqc .
```

Compare pre- and post-trim quality to confirm adapter removal and tail improvement.

---

## 🚀 Next Steps

You can now proceed to:

1. **Quantification** using *Salmon* (alignment-free)
   or
2. **Alignment + Counting** using *STAR* and *featureCounts*.

---

**Author:** *Imran Hasan*
**Cluster:** *Negishi (Purdue University)*
**Date:** *November 2025*

```

---


```

