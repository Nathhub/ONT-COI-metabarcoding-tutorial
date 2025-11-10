# COI Metabarcoding Data Processing with Oxford Nanopore PromethION  
### OceanX Around Africa Expedition  
*A step-by-step teaching guide for generating OTU and taxonomy tables from COI amplicons.*

---

## 🧭 Introduction

This tutorial guides you through the full bioinformatics workflow used to process **COI metabarcoding data** generated with the **Oxford Nanopore PromethION** platform during the *OceanX Around Africa* expedition.  

The goal is to produce two main outputs:
1. An **OTU table** — showing read abundances of each consensus sequence (OTU) across samples.  
2. A **taxonomy table** — linking each OTU to its best taxonomic match from a reference database.

Each step below explains **what**, **why**, and **how**, so that students can both *replicate* and *understand* the reasoning behind each command.

---

## 🧩 Step 1 — Combine FASTQ Files by Barcode

**Purpose:**  
The Nanopore sequencer produces multiple FASTQ files per barcode (one per read batch). We merge them to have one file per barcode for easier downstream processing.

**Command:**
```bash
seqkit scat -j 32 -f folder_with_fastq/ > combined.fastq.gz
```

**Notes:**  
- `seqkit scat` merges all FASTQ files in the folder.  
- `-j 32` uses 32 threads for faster execution.  
- Each `combined.fastq.gz` file corresponds to one barcode (sample).

---

## 🔍 Step 2 — Assess Read Quality and Length

**Purpose:**  
Before trimming and filtering, it’s good practice to inspect read quality and length distribution using **NanoPlot**.

**Command:**
```bash
NanoPlot --fastq combined.fastq.gz --threads 32 --outdir nanoplot_reports/
```

**Notes:**  
- Look for the average read length and quality scores (Q-scores).  
- For COI amplicons, you expect lengths around 300–350 bp.

---

## ✂️ Step 3 — Primer Trimming with Cutadapt

**Purpose:**  
Remove primer sequences from both ends of each read. This ensures that downstream clustering and alignment are based on the biological insert only.

**Command (for COI):**
```bash
cutadapt    -g GGWACWGGWTGAACWGTWTAYCCYCC...TGRTTYTTYGGICAYCCIGARGTITA  \
            -g TAIACYTCIGGRTGICCRAARAAYCA...GGRGGRTAWACWGTTCAWCCWGTWCC   \
            --error-rate 0.2   \
            --cores=64  \
            --length 350   \
            --discard-untrimmed  \
            -o trimmed_barcode${i}.fastq.gz  \
            barcode${i}.fastq.gz
```

**Notes:**  
- The two `-g` options represent both primer orientations (forward and reverse complement).  
- `--error-rate 0.2` allows up to 20% mismatch due to ONT error rates.  
- `--discard-untrimmed` removes reads without recognizable primer pairs.

---

## 🧹 Step 4 — Quality and Length Filtering with Chopper

**Purpose:**  
Filter out reads that are too short, too long, or have low mean quality scores.

**Command:**
```bash
gunzip -c trimmed_barcode${i}.fastq.gz | chopper     -q 22     --minlength 290 --maxlength 350     -t 64 | gzip > filtered_barcode${i}.fastq.gz
```

**Notes:**  
- `-q 22` keeps reads with mean quality ≥ Q22.  
- Length thresholds (290–350 bp) correspond to the COI fragment size.

---

## 🧬 Step 5 — Generate Consensus Sequences (Amplicon Sorter)

**Purpose:**  
Cluster reads from each sample into groups of identical or near-identical sequences and generate consensus sequences (representing potential OTUs).

**Command:**
```bash
python3 amplicon_sorter/amplicon_sorter.py     -min 290     -i filtered_barcode${i}.fastq.gz     -o sorted_barcode${i}     -np 64 -ra -maxr 1000000
```

**Notes:**  
- `-ra` enables random subsampling (useful for large datasets).  
- The output folder contains `consensusfile.fasta` for each barcode.

---

## 🧱 Step 6 — Combine and Deduplicate Consensus Sequences

**Purpose:**  
Merge all per-sample consensus sequences into one FASTA file and remove identical duplicates.

**Commands:**
```bash
cat sorter_out/sorted_barcode{01..18}/consensusfile.fasta > all_samples_consensus.fasta
```

**Remove identical sequences:**
```bash
awk '
BEGIN {seq=""; name=""}
(/^>/) {
    if (seq != "") {
        if (!(seq in seen)) {
            seen[seq] = name; order[++i] = seq
        }
        seq = ""
    }
    name = $0
}
!/^>/ { seq = seq $0 }
END {
    if (!(seq in seen)) {
        seen[seq] = name; order[++i] = seq
    }
    for (j=1; j<=i; j++) {
        print seen[order[j]]; print order[j]
    }
}' all_samples_consensus.fasta > all_samples_consensus_unique.fasta
```

**Remove duplicate headers:**
```bash
awk '/^>/{name=$0; if(seen[name]++) next} {print}' all_samples_consensus_unique.fasta > all_samples_consensus_nodup.fasta
```

---

## 🧾 Step 7 — Prepare Reads for Mapping

**Purpose:**  
Convert filtered reads to FASTA and standardize headers so that each read has a unique ID containing the sample name.

**Commands:**
```bash
seqkit fq2fa filtered_reads.fastq.gz -o filtered_reads.fasta
```

**Rename headers:**
```bash
for i in {01..18}; do
  sample="barcode$i"
  input="filtered_${sample}.fasta"
  output="renamed_${sample}.fasta"

  awk -v sample="$sample" '
    BEGIN {i=0}
    /^>/ {
      i++
      header = sample "_" i
      print ">" header ";sample=" sample
      next
    }
    { print }
  ' "$input" > "$output"
done
```

---

## 📊 Step 8 — Construct the OTU Table with Minimap2

**Purpose:**  
Map all reads back to the consensus sequences (OTUs) to obtain per-sample read counts.

### 1️⃣ Index the reference (consensus sequences)
```bash
minimap2 -d consensus.mmi all_samples_consensus_nodup.fasta
```

### 2️⃣ Map reads from each sample
```bash
for f in fasta_reads/*.fasta; do
  base=$(basename "$f" .fasta)
  minimap2 -ax map-ont consensus.mmi "$f" > "${base}.sam"
done
```

### 3️⃣ Convert SAM → BAM, sort, and index
```bash
for f in *.sam; do
  base=$(basename "$f" .sam)
  samtools view -bS "$f" | samtools sort -o "${base}.sorted.bam"
  samtools index "${base}.sorted.bam"
done
```

### 4️⃣ Keep primary alignments only
```bash
mkdir -p primary_bams
parallel -j 30 --will-cite '
  b={};
  base=$(basename "$b" .sorted.bam);
  samtools view -b -F 0x904 "$b"     | samtools sort -o primary_bams/${base}.primary.sorted.bam
  samtools index primary_bams/${base}.primary.sorted.bam
' ::: *.sorted.bam
```

### 5️⃣ Generate per-sample counts
```bash
samples=($(ls -1 primary_bams/*.primary.sorted.bam | sed 's#.*/##; s/\.primary\.sorted\.bam$//' | sort))

mkdir -p counts
parallel -j 30 --will-cite '
  f={};
  base=$(basename "$f" .primary.sorted.bam);
  samtools idxstats "$f"     | awk -v OFS="\t" '"'"'$1!="*" {print $1,$3}'"'"'     > counts/"$base".tsv
' ::: primary_bams/*.primary.sorted.bam
```

### 6️⃣ Assemble the OTU table
```bash
printf "OTU" > otu_table.primary.tsv
for s in "${samples[@]}"; do printf "\t%s" "$s" >> otu_table.primary.tsv; done
printf "\n" >> otu_table.primary.tsv

grep '^>' all_samples_consensus_nodup.fasta | sed 's/^>//' > otu_order.txt

awk -v OFS="\t" -v SAMPLES="$(printf "%s " "${samples[@]}")" '
  BEGIN{
    n=split(SAMPLES, samp, " ");
    for(i=1;i<=n;i++){
      file="counts/" samp[i] ".tsv";
      while((getline < file)>0){ C[$1, samp[i]]=$2+0; seen[$1]=1 }
      close(file);
    }
  }
  FNR==NR { order[++m]=$0; next }
  END{
    for(i=1;i<=m;i++){
      otu=order[i];
      printf "%s", otu;
      for(j=1;j<=n;j++){
        printf "\t%s", ((otu SUBSEP samp[j]) in C ? C[otu, samp[j]] : 0);
      }
      printf "\n";
    }
  }
' otu_order.txt >> otu_table.primary.tsv
```

---

## ⚖️ Step 9 — Normalize Counts (CPM)

**Purpose:**  
Normalize counts to **Counts Per Million (CPM)** to account for differences in sequencing depth across samples.

**Command:**
```bash
awk 'NR==1{print; next}
     {sum=0; for(i=2;i<=NF;i++) sum+=$i;
      printf "%s", $1;
      for(i=2;i<=NF;i++) printf "\t%.6f", (sum?($i/sum*1e6):0);
      printf "\n"}' otu_table.primary.tsv > otu_table.primary.cpm.tsv
```

---

## 🧭 Step 10 — Assign Taxonomy Using BLAST

**Purpose:**  
Identify each OTU by comparing it to a reference COI database (here, **MZGdb_COI**).

**Commands:**
```bash
makeblastdb -in ~/ote_db/MZG_db/all/MZGfasta-coi__MZGdbALL__o00__A.fasta             -dbtype nucl             -out MZGdb_COI

blastn -query all_samples_consensus_nodup.fasta        -db MZGdb_COI        -out blast_consensus_nodup_results.b6        -outfmt 6        -evalue 1e-10        -num_threads 64        -max_target_seqs 1
```

**Notes:**  
- The output format (`-outfmt 6`) is tabular and easy to parse.  
- `-max_target_seqs 1` ensures you keep only the top hit per OTU.  

---

## 🌳 Step 11 — Build a Taxonomy Table

**Purpose:**  
Combine BLAST results with taxonomy metadata and format hierarchical taxonomic levels.

### Join BLAST hits to taxonomy strings:
```bash
awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR { tax[$1]=$2; next }
{
  t = ($2 in tax ? tax[$2] : "")
  print $1, $2, t
}' ~/ote_db/MZG_db/all/MZGmothur-coi__MZGdbALL__o00__A.txt blast_consensus_nodup_results.b6 | sort -k1,1 > best_hits_with_tax.tsv
```

### Determine number of taxonomic levels:
```bash
MAXLVL=$(awk -F'\t' '{
  gsub(/;$/,"",$3);
  n=split($3,a,/;/);
  if(n>max) max=n
} END{print (max?max:0)}' best_hits_with_tax.tsv)
```

### Create header and expand taxonomy columns:
```bash
{
  printf "QueryID\tSubjectID"
  for ((i=1; i<=MAXLVL; i++)); do
    printf "\tL%d" "$i"
  done
  printf "\n"
} > taxonomy_table.tsv

awk -v MAXLVL="$MAXLVL" -F'\t' 'BEGIN{OFS="\t"}
{
  gsub(/;$/,"", $3)
  n=split($3, a, /;/)
  printf "%s\t%s", $1, $2
  for(i=1;i<=MAXLVL;i++){
    printf "\t%s", (i<=n ? a[i] : "")
  }
  printf "\n"
}' best_hits_with_tax.tsv >> taxonomy_table.tsv
```

### Rename columns for readability:
```text
"Kingdom", "Subkingdom_1", "Subkingdom_2", "Phylum", "Subphylum", "Subphylum_1", 
"Superclass", "Class", "Subclass", "Infraclass", "Superorder", "Order", 
"Suborder_1", "Suborder_2", "Suborder_3", "Family", "Subfamily", 
"Genus", "Subgenus", "Species"
```

---

## 📘 References

- **Cutadapt** – Martin, M. (2011). *EMBnet.journal*, 17(1):10–12.  
- **Chopper** – Oxford Nanopore read quality and length filter.  
- **Amplicon Sorter** – clustering and consensus generation for amplicon reads.  
- **Minimap2** – Li, H. (2018). *Bioinformatics*, 34(18):3094–3100.  
- **Samtools** – Danecek et al. (2021). *Gigascience*, 10(2).  
- **MZG COI Database** – curated COI reference sequences for metazoan identification.

---

**End of tutorial.**  
🧬 *You now have a complete OTU table (raw and normalized) and a taxonomy table ready for ecological or diversity analyses.*
