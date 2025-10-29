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











module --force purge

module load biocontainers

module load fastqc/0.11.9

vdb-config --prefetch-to-cwd

module load fastqc
