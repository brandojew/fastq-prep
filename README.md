# fastq-prep

Prepares sequencing data for parallelized variant calling. Takes BAM or FASTQ files with paired-end reads and produces split and compressed FASTQ files. Corresponding reads occur on same line in respective FASTQ file. 


Currently working on allowing many types of input (BAM, paired and interleaved FASTQ, compressed files, etc) and the option for producing paired end or interleaved chunks. Also, filtering of reads (secondary, qc-fail, etc) is not included. 

Usage:
  * python fastq_prep.py [output.prefix] [input.bam or input.fq]
 OR python fastq_prep.py [output.prefix] [input_R1.fq,input_R2.fq]

  * Multiple input files are separated by commas, no spaces

Requires:
  * pysam-0.7.5 (Really old version, working on compatibility with current verisons)


