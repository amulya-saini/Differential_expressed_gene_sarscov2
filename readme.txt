Assignment 2
Date: 12/03/2023
Subject: B636: Genomic Data Analytics and Precision Medicine


Unix Commands used:

###### download SRA files ###### SRA TOOLKIT

- load the module
module load sra-toolkit

-Download the data
prefetch ascension_numbers
fastq-dump ascension_numbers

###### Quality control ###### FASTQC , MULTIQC

- loading fastqc
module load fastqc
fastqc path_to_the_fastq_files

- combining all the fastqc reports using multiqc
- Install the module
module load python
pip install multiqc

- add the module to path (in .bashrc) 
- navigate to the folder in which all the fastqc reports are present using cd
- once in the folder run 
multiqc .

###### Trimming ###### TRIMGALORE

- loading module
module load python
module load trimgalore

- Trimming using trim_galore with said parameters (Phred score 25, minimum length 25)
trim_galore --quality 25 --length 16 --stringency 3 --illumina --fastqc --output_dir trimmed_data *.fastq

- running the fastqc and multiqc 

###### Alignment ###### BOWTIE2

- loading module
module load bowtie2

- indexing of reference genome
bowtie2-build reference_genome.fa bowtie2_index/index

- alignment of every sample using bowtie2 index

bowtie2_index="path/to/bowtie2_index"

# Specify the location of the trimmed fastq files
trimmed_data_dir="path/to/trimmed_data"

# Specify the location to store SAM files
sam_files_dir="path/to/store/samfiles"

# Iterate over the range of file numbers SRR22269872 - 83
for i in {72..83}; do
    # Construct the paths for input and output files
    input_fastq="${trimmed_data_dir}/SRR222698${i}_trimmed.fq"
    output_sam="${sam_files_dir}/SRR222698${i}.sam"

    # Run Bowtie2 alignment
    bowtie2 -x "${bowtie2_index}" -U "${input_fastq}" -S "${output_sam}"
done

###### Sam to Bam conversion ###### SAMTOOLS

- loading modules
module swap gcc/12.1.0 gcc/9.3.0
module load python
module load samtools

# Specify the directory where SAM files are located
sam_dir="path/to/the/sam_files"

# Specify the directory to store BAM files
bam_dir="path/to/store/bam_files"

# Navigate to the directory containing SAM files
cd "$sam_dir" || exit

# Iterate over the range of file numbers SRR22269872 - 83
for i in {72..83}; do
    sam_file="SRR222698${i}.sam"

    # Generate the output BAM file name in the BAM directory
    bam_file="$bam_dir/SRR222698${i}.bam"

    # Converting SAM to BAM
    samtools view -bS -o "$bam_file" "$sam_file"

    # sorting and indexing the BAM files
    samtools sort -o "${bam_file%.bam}_sorted.bam" "$bam_file"
    samtools index "${bam_file%.bam}_sorted.bam"
done

###### Quantification ###### FEATURECOUNTS

- loading the module
module load subread

# Specify the path to the folder containing BAM files
bam_folder="path/to/the/folder/containing/sorted_bam_files"

# Specify the paths to gene and miRNA annotation files
gene_gtf_file="/path/to/the/gene_annotation.gtf"
miRNA_gff_file="/path/to/the/miRNA_annotation.gff3"

# Specify the output folder
output_folder="/path/to/store/the/counts"

# Run featureCounts for gene annotations
gene_output_file="$output_folder/all_samples_gene_count.txt"
featureCounts -a "$gene_gtf_file" -t gene -g 'gene_id' -o "$gene_output_file" "${bam_folder}"/*.bam

# Run featureCounts for miRNA annotations
miRNA_output_file="$output_folder/all_samples_miRNA_count.txt"
featureCounts -a "$miRNA_gff_file" -t miRNA,miRNA_primary_transcript -g 'Name' -o "$miRNA_output_file" "${bam_folder}"/*.bam

##################Moving the files to the local directory for further downstream analysis in R (DEseq2) #######################
