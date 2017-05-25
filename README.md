# fastq-prep

Prepares sequencing data for parallelized variant calling. Takes BAM, SAM, CRAM, or FASTQ files with paired-end reads and produces split and compressed FASTQ files. Corresponding reads occur on same line in respective FASTQ file. 


Currently filters reads with SAM flags indicating non primary alignment and quality control failure. Also filters reads with SAM flags indicating "first in pair" or "second in pair" when such a read has already been encountered (i.e. Two reads have the same name and are both reported as first in pair). Filtering options are specified in the header of the python file. 

### usage
  * python fastq_prep.py [output.prefix] [input.bam or input.fq] 
  * OR python fastq_prep.py [output.prefix] [input_R1.fq,input_R2.fq]
  * OR [write SAM records to stdout] | python fastq_prep.py [output.prefix]

  * Multiple input files are separated by commas, no spaces

### requirements
  * pysam

### todo
  * Support generation of interleaved FASTQ chunks.
