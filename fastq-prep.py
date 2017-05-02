import pysam
import gzip
import sys


#Global constants for development
RECORDS_PER_FILE = 750000 # ~60MB
COMPRESS_LEVEL = 5
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

class RecordWriter(object):
  """Interface for writing records to split paired-end FASTQ files
    
    Formats records to FASTQ format and writes to files in chunks of 
    roughly 60 MB for parallelization of downstream analysis. Unpaired
    records are written to one file.

  Attributes:
    fastq_prefix: Full path prefix for FASTQ files generated. Index string
      and file extension will be appended to this.
    fastq_index: Current index of chunk being written to.
    fastq_records: Number of records written to current chunk.
    fastq_file_1: Current chunk for first-in-pair records.
    fastq_file_2: Current chunk for second-in-pair records.
    unpaired_file: File to write all unpaired records to.
  """
  def __init__(self, fastq_prefix):
    """Inits RecordWriter for given FASTQ prefix"""
    self.fastq_prefix = fastq_prefix
    self.fastq_index = 0
    self.fastq_records = 0
    self.fastq_file_1 = None
    self.fastq_file_2 = None
    self.update_fastq_index()
    unpaired_path = ''.join([fastq_prefix,'_unpaired.fastq.gz'])
    self.unpaired_file = gzip.open(unpaired_path, 'w', compresslevel = 5)

  def update_fastq_index(self):
    """Closes current chunks for FASTQ files and creates next chunks"""
    if self.fastq_file_1 and self.fastq_file_2:
      self.fastq_file_1.close()
      self.fastq_file_2.close()

    self.fastq_index += 1
    self.fastq_records = 0
    index_str = str(self.fastq_index).zfill(5)

    fastq_path_1 = \
      '{}_{}_{}.fastq.gz'.format(self.fastq_prefix, index_str, "R1")
    fastq_path_2 = \
      '{}_{}_{}.fastq.gz'.format(self.fastq_prefix, index_str, "R2")
    self.fastq_file_1 = gzip.open(fastq_path_1, 'w', compresslevel = COMPRESS_LEVEL)
    self.fastq_file_2 = gzip.open(fastq_path_2, 'w', compresslevel = COMPRESS_LEVEL)

  def write_paired_records(self, record_a, record_b):
    """Writes given records to appropriate FASTQ chunk
    
    Determines which read is 1st in pair, and whether each sequence maps to
    reverse strand. Writes appropriately formatted records to respective chunk

    Args:
      record_a, record_b: Paired records (pysam.Samfile.fetch() item)
      in no particular order
    """
    if record_a.is_read1:
      record_1 = record_a
      record_2 = record_b
    else:
      record_1 = record_b
      record_2 = record_a

    if record_1.is_reverse:
      self.fastq_file_1.write("\n".join([record_1.qname+"/1", \
        reverse_complement(record_1.seq),"+", \
        record_1.qual])+"\n")
    else:
      self.fastq_file_1.write\
        ("\n".join([record_1.qname+"/1",record_1.seq,"+",record_1.qual])+"\n")

    if record_2.is_reverse:
      self.fastq_file_2.write("\n".join([record_2.qname+"/2", \
        reverse_complement(record_2.seq),"+", \
        record_2.qual])+"\n")
    else:
      self.fastq_file_2.write\
        ("\n".join([record_2.qname+"/2",record_2.seq,"+",record_2.qual])+"\n")
    self.fastq_records += 1

    if self.fastq_records == RECORDS_PER_FILE:
      self.update_fastq_index()

  def write_unpaired_record(self, record):
    """Writes given record to unpaired FASTQ file and formats appropriately"""
    if record.is_reverse:
      self.unpaired_file.write("\n".join(\
        [record.qname,reverse_complement(record.seq),"+",record.qual])+"\n")
    else:
      self.unpaired_file.write\
        ("\n".join([record.qname,record.seq,"+",record.qual])+"\n")

class SimpleRecord(object):
  """Object that stores the record entries needed for a FASTQ record
  
  Attributes:
    qname: Name of record
    seq: Sequence string
    qual: Quality string
  """
  def __init__(self, qname, seq, qual):
    """Inits object with FASTQ entries"""
    self.qname = qname
    self.seq = seq
    self.qual = qual

def reverse_complement(seq):
  """Returns reverse complement of given sequence using list comprehension"""
  return "".join([COMPLEMENT[base] for base in seq[::-1]]) 

def bam_to_fastq(bam_path, fastq_prefix):
  """Converts BAM file into split FASTQ files

  Args:
    bam_path: Full path to BAM file to convert.
    fastq_prefix: Full path prefix for FASTQ files to generate.
  """
  bam_file = pysam.Samfile(bam_path, "rb")
  record_writer = RecordWriter(fastq_prefix)
  record_count = 0
  record_dict = {}
  # TODO: If memory is an issue, consider pickling long_term_record_dict until 
  #       unmapped reads reached (Check if BAM sorted by coordinate). 
  """long_term_record_dict = {}"""

  for record in bam_file.fetch():
    record_count += 1
    if record.is_paired:
      if record.qname in record_dict:
        record_writer.write_paired_records(record,record_dict.pop(record.qname))
        # Implementation of long_term_record_dict population
        """
      elif not record.is_unmapped and record.mate_is_unmapped:
        #Record is mapped but mate is not, store in long-term
        long_term_record_dict[record.qname] = record
      elif record.is_unmapped and not record.mate_is_unmapped:
        #Record is unmapped but mate is mapped, so mate stored in long_term
        if record.qname in long_term_record_dict:
          record_writer.write_paired_records\
            (record,long_term_record_dict.pop(record.qname))
        else:
          #Mate does not exist
          record_writer.write_unpaired_record(record)
        """
      else:
        record_dict[record.qname] = record
    else:
      record_writer.write_unpaired_record(record)
    if not record_count % 10000000:
      print('Processed {} records'.format(str(record_count)))

  for record in record_dict.values():
    record_writer.write_unpaired_record(record)


def read_fastq_record(fastq_file):
  """Reads one record from a FASTQ file

  Args:
    fastq_file: Open file object
  Returns:
    record: Object with qual, seq, qname attributes
  """
  qname = fastq_file.readline().strip()
  if not qname or qname[0] != "@":
    return
  seq = fastq_file.readline().strip()
  if fastq_file.readline()[0] != "+":
    return 
  qual = fastq_file.readline().strip()
  return SimpleRecord(qname, seq, qual)


def split_existing_interleaved_fastq():
  try:
    input_fastq = gzip.open(fastq_path, 'r')
    input_fastq.readline()
  except IOError:
    input_fastq = open(fastq_path, 'r')
  input_fastq.seek(0)
  # TODO: Read file and submit adjacent records to RecordWriter

def split_existing_paired_fastq():
  try:
    input_fastq_1 = gzip.open(fastq_1_path, 'r')
    input_fastq_2 = gzip.open(fastq_1_path, 'r')
    input_fastq_1.readline()
    input_fastq_2.readline()
  except IOError:
    input_fastq_1 = open(fastq_1_path, 'r')
    input_fastq_2 = open(fastq_2_path, 'r')
  input_fastq_1.seek(0)
  input_fastq_2.seel(0)
  # TODO: Read files and submit parallel records to RecordWriter

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python bam_to_fastq.py [in.bam] [out.fastq prefix]")
    sys.exit()
    
  bam_path = sys.argv[1]
  fastq_prefix = sys.argv[2]
  # TODO: Check input for which function to call

  print('Converting following BAM to FASTQ format:')
  print(bam_path)
  print('FASTQ files will be split into files with the following structure:')
  print('{}_XXXXX_RX.fastq.gz'.format(fastq_prefix))

  bam_to_fastq(bam_path, fastq_prefix)

  print("Done")
