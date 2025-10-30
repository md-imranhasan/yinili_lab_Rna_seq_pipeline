# yinili_lab_Rna_seq_pipeline

rna-seq/
├─ README.md
├─ envs/
│  └─ rnaseq.yml
├─ config/
│  ├─ samples.csv            # sample metadata
│  └─ design.yml             # contrasts, organism, strandedness
├─ data/
│  ├─ fastq/                 # raw .fastq.gz
│  └─ reference/             # genome, GTF, Salmon index
├─ results/
│  ├─ qc/fastqc/
│  ├─ qc/multiqc/
│  ├─ trim/
│  ├─ align/                 # STAR/HISAT2 BAMs (if alignment path)
│  ├─ quantify/              # Salmon/Kallisto outputs (if lightweight path)
│  ├─ counts/                # featureCounts matrices
│  ├─ de/                    # DESeq2 outputs
│  └─ enrich/                # GO/KEGG/GSEA outputs
└─ workflow/
   ├─ Snakefile              # optional
   └─ scripts/
      ├─ make_tx2gene.R
      └─ deseq2.R




Perfect 👍 Let’s slow it down and make it **super clear and simple**.
You already have:

* A **metadata file** (Accession list) for *Frontal Cortex – ALS spectrum MND*
* A **metadata file** (Accession list) for *Frontal Cortex – non-neurological control*
  and you’re working on a cluster using the **SRA-Toolkit module**.

---

## 🧩 Step-by-Step Explanation

### 1️⃣ Create and move to your project folder

You want one main folder for this dataset:

```bash
cd /depot/yinili/data/Li_lab/
mkdir -p GSE124439_Hammell2019/Frontal_Cortex
cd GSE124439_Hammell2019/Frontal_Cortex
```

Now you’re inside your working folder:

```
/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex
```

---

### 2️⃣ Load SRA-Toolkit

You must load the right software first (this adds the commands like `prefetch`, `fasterq-dump`, etc. to your path):

```bash
module --force purge
module load biocontainers
module load sra-tools/2.11.0-pl5262
```

---

### 3️⃣ Configure where downloaded data will go

Run this once to make sure SRA-Toolkit downloads files *in your current folder*:

```bash
vdb-config --prefetch-to-cwd
```

That’s it — it will save `.sra` files right here.

---

### 4️⃣ Place your Accession list here

Put your text files (each line = one SRR ID) in this folder.
Example:

**File 1:** `ALS_list.txt`

```
SRR11941281
SRR11941282
SRR11941283
```

**File 2:** `CTRL_list.txt`

```
SRR11941284
SRR11941285
```

So now you have:

```
/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/ALS_list.txt
/depot/yinili/data/Li_lab/GSE124439_Hammell2019/Frontal_Cortex/CTRL_list.txt
```

---










module --force purge

module load biocontainers
module load sra-tools/2.11.0-pl5262

module load fastqc/0.11.9

vdb-config --prefetch-to-cwd

module load fastqc
module load multiqc






#SRA fastq download

while read SRR; do echo "Downloading $SRR ..."; fasterq-dump "$SRR" -O fastq/; done < ALS_list.txt

I want to use threads 32

