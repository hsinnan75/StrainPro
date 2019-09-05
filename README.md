StrainPro -- a highly accurate Metagenomic strain-level profiling tool
==============================
Developers: Dr. Hsin-Nan Lin, Dr. Yaw-Ling Lin and Dr. Wen-Lian Hsu Institute of Information Science, Academia Sinica, Taiwan.

# Introduction
Characterizing the taxonomic diversity of a microbial community is very important to understand the roles of microorganisms. Next generation sequencing (NGS) provides great potential for investigation of a microbial community and leads to Metagenomic studies. NGS generates DNA sequences directly from microorganism samples, and it requires analysis tools to identify microbial species (or taxonomic composition) and estimate their relative abundance in the studied community. Here we developed a novel metagenomic analysis tool, called StrainPro, which is highly accurate both at characterizing microorganisms at strain-level and estimating their relative abundances. A unique feature of StrainPro is it identifies representative sequence segments from reference genomes. We generate three simulated datasets using known strain sequences and another three simulated datasets using unknown strain sequences.

# Download
Please use the command
  ```
  $ git clone https://github.com/hsinnan75/StrainPro.git
  ```
to download the package of StrainPro.

# Dependencies
To compile StrainPro, it requires libboost-all-dev, libbz2-dev, and liblzma-dev installed in your system.

# Compiling
Please change to StrainPro's folder and type 'make' to compile all programs of StrainPro.

# Taxonomy & Reference genomes
To download taxonomy information, please use the script "download_taxonomy.sh" to download taxonomy files. (requisite for StrainPro)
  ```
  $./download_taxonomy.sh
  ```

To download reference genomes, please use the script "download_genomic_library.sh" to download reference genomes.
  ```
  $./download_genomic_library.sh library (library: archaea | bacteria | viral | fungi | human)
  ```
The script will download all complete chromosomes/genomes of the designated library.

If you would like to use a customized reference library, please make sure your reference genomes follow the below requirements:
 - Genome sequences are in FASTA format.
 - Each sequence's header must start with "taxid|xxxxx|" where xxxxx is the taxon ID. For example,
   >taxid|562|NZ_CP027599.1 Escherichia coli strain 97-3250 
   ATCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGT....

# Database
To build a database, please use the StrainPro-build command
  ```
  $./StrainPro-build -r reference-fna -o ref_idx [ref_idx is the output folder for BWT indexes].
  ```
It will take a while to process the genome sequences and build BWT indexes.

# Mapping & Profiling
To characterize the taxonomic composition of the input metagenomic data, please use the StrainPro-map command
  ```
  $./StrainPro-map -i idx -f read -o profile.txt [idx=BWT_index prefix / folder of BWT_indexes, read=NGS_dataset]
  ```
You may specify a specific index or a folder of multiple indexes. A specific index is designated with its prefix. The input NGS data is in FASTQ or FASTA format. Multiple files are separated with a space character. The output is a text file.

# Output format
StrainPro ouputs the taxonomic composition of the input metagenomic data directly. The output format is
  ```
TaxID  Read_count  Est_depth   Est_relative_abundance   Confidence_score
  ```
where TaxID is the NCBI taxon identifier; Read_count is the number of reads that are classified into that taxon id; Est_depth is the estimated read depth of that taxon id; Est_relative_abundance is the estimated relative abundance (percentage); Confidence_score is the confidence score of that prediction.

