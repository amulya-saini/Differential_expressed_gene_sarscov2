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

I used DEseq2 to perform differentially expressed genes in the attached R code. 
- Drop columns ("Chr", "Start", "End", "Strand", "Length") from the counts' files obtained from feature counts
- Change the first column to row names
- Combine both the datasets without losing row/column values
- filter the data frame by summing the counts across each and dropping rows(genes) with less than 1 count.
- Create DEseq dataset from the matrix using the meta data from the study
- Run DEseq Analysis
- Sort the adjusted p-values in ascending order
- Gene with less than 0.05 adjusted p_value are selected as Differentially expressed genes


Data used:

- I have used Human reference genome version 38 (hg38), both fasta file and gtf file were downloaded from the following website(illumina):
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz

- I have downloaded the miRNA annotation gff3 file from the following website (miRBase):
https://mirbase.org/download/hsa.gff3

- The RNA seq data was obtained from the following study
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA901149


