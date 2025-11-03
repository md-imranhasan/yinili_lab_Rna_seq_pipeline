Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: RNA-Seq Pipeline


Perfect 💪 — below is a **ready-to-upload, documentation-style Markdown (`README.md`)** that captures your full RNA-seq workflow **from raw SRA download to trimming with Trimmomatic** on the Negishi HPC cluster.

You can copy–paste this into a file named `README.md` in your project or GitHub repo.

---

```markdown
# 🧬 RNA-seq Processing Pipeline — GSE124439 (Hammell 2019)

This repository documents a full end-to-end workflow for processing the **GSE124439** RNA-seq dataset  
(`PRJNA512012`) using the Negishi HPC cluster at Purdue.  
The example focuses on the **Frontal Cortex** and **Motor Cortex (Medial)** sub-regions.

---

## 📁 Directory Layout

```

/depot/yinili/data/Li_lab/GSE124439_Hammell2019/
├── adapters/
│   └── TruSeq3-PE.fa
├── Frontal_Cortex/
│   ├── case/
│   │   ├── fastq/
│   │   ├── qc/
│   │   ├── trim_trimmomatic/
│   │   └── logs/
│   ├── control/
│   └── SraRunTable_Frontal_Cortex.csv
└── motor_cortex_(medial)/
└── case/

````

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

---

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

**Author:** *Asif Ahmed*
**Cluster:** *Negishi (Purdue University)*
**Date:** *November 2025*

```

---

Would you like me to extend this same file with the **next stage (Salmon quantification)** section so you can keep one continuous GitHub documentation?
```

