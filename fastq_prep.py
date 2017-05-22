"""Converts BAM or FASTQ files to split and compressed FASTQ files"""

import gzip
import sys
import os
import pysam

#Global constants for development
RECORDS_PER_FILE = 750000 # ~60MB
COMPRESS_LVL = 5
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
    unpaired_path = ''.join([fastq_prefix, '_unpaired.fastq.gz'])
    self.unpaired_file = gzip.open(unpaired_path, 'wt',
                                   compresslevel=COMPRESS_LVL)

  def update_fastq_index(self):
    """Closes current chunks for FASTQ files and creates next chunks"""
    if self.fastq_file_1 and self.fastq_file_2:
      self.fastq_file_1.close()
      self.fastq_file_2.close()
    self.fastq_index += 1
    self.fastq_records = 0
    index_str = str(self.fastq_index).zfill(5)
    fq_path_1 = '{}_{}_{}.fastq.gz'.format(self.fastq_prefix, index_str, "R1")
    fq_path_2 = '{}_{}_{}.fastq.gz'.format(self.fastq_prefix, index_str, "R2")
    self.fastq_file_1 = gzip.open(fq_path_1, 'wt', compresslevel=COMPRESS_LVL)
    self.fastq_file_2 = gzip.open(fq_path_2, 'wt', compresslevel=COMPRESS_LVL)

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
    qname_1 = '@{}/1'.format(record_1.qname)
    qname_2 = '@{}/2'.format(record_2.qname)
    if record_1.is_reverse:
      self.fastq_file_1.write("\n".join([qname_1,
                                         reverse_complement(record_1.seq),
                                         "+", record_1.qual, ""]))
    else:
      self.fastq_file_1.write("\n".join([qname_1,
                                         record_1.seq, "+",
                                         record_1.qual, ""]))
    if record_2.is_reverse:
      self.fastq_file_2.write("\n".join([qname_2,
                                         reverse_complement(record_2.seq),
                                         "+", record_2.qual, ""]))
    else:
      self.fastq_file_2.write("\n".join([qname_2,
                                         record_2.seq, "+",
                                         record_2.qual, ""]))
    self.fastq_records += 1
    if self.fastq_records == RECORDS_PER_FILE:
      self.update_fastq_index()

  def write_unpaired_record(self, record):
    """Writes given record to unpaired FASTQ file and formats appropriately"""
    qname = '@{}'.format(record.qname)
    if record.is_reverse:
      self.unpaired_file.write("\n".join(
        [qname, reverse_complement(record.seq), "+", record.qual, ""]))
    else:
      self.unpaired_file.write("\n".join(
        [qname, record.seq, "+", record.qual, ""]))

class SimpleRecord(object):
  """Object that stores the record entries needed for a FASTQ record

  Attributes:
    qname: Name of record
    seq: Sequence string
    qual: Quality string
  """
  def __init__(self, qname, seq, qual, is_read1):
    """Inits object with FASTQ entries"""
    self.qname = qname
    self.seq = seq
    self.qual = qual
    self.is_read1 = is_read1
    self.is_reverse = False

  def copy_constructor(self, record):
    """Inits object using Pysam record object"""
    self.qname = record.qname
    self.seq = record.seq
    self.qual = record.qual
    self.is_read1 = record.is_read1
    self.is_reverse = record.is_reverse

  def fastq_format(self):
    """Prints data members in FASTQ format"""
    print("\n".join([self.qname, self.seq, "+", self.qual, ""]))

def reverse_complement(seq):
  """Returns reverse complement of given sequence"""
  return "".join([COMPLEMENT[base] for base in seq[::-1]])

def split_bam(bam_path, fastq_prefix):
  """Converts BAM or SAM file into split FASTQ files

  Args:
    bam_path: Full path to BAM or SAM file to convert.
    fastq_prefix: Full path prefix for FASTQ files to generate.
  """
  try:
    bam_file = pysam.Samfile(bam_path, "rb")
  except ValueError:
    bam_file = pysam.Samfile(bam_path, "r")
  record_writer = RecordWriter(fastq_prefix)
  record_count = 0
  record_dict = {}
  for record in bam_file:
    record_count += 1
    if record.is_paired:
      if record.qname in record_dict:
        record_writer.write_paired_records(
          record, record_dict.pop(record.qname))
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
    record: Object with qual, seq, qname attributes.
      Will return None when EOF reached (or empty line)
  """
  qname = fastq_file.readline().strip()
  if not qname:
    return
  if qname[0] != "@":
    raise Exception("Invalid FASTQ entry")
  else:
    qname = qname[1::]
  if qname[-1] == "1":
    is_read1 = True
  elif qname[-1] == "2":
    is_read1 = False
  else:
    is_read1 = None
  seq = fastq_file.readline().strip()
  if fastq_file.readline()[0] != "+":
    raise Exception("Invalid FASTQ entry")
  qual = fastq_file.readline().strip()
  if (qname[len(qname)-2:len(qname)] == '/1'
      or qname[len(qname)-2:len(qname)] == '/2'):
    qname = qname[0:len(qname)-2]
  return SimpleRecord(qname, seq, qual, is_read1)

def split_interleaved_fastq(fastq_path, fastq_prefix):
  """Converts a single interleaved FASTQ file to split and paired chunks

  Args:
    fastq_path: Full path to input interleaved FASTQ file.
    fastq_prefix: Prefix for output FASTQ chunks
  """
  try:
    input_fastq = gzip.open(fastq_path, 'r')
    input_fastq.readline()
  except IOError:
    input_fastq = open(fastq_path, 'r')
  input_fastq.seek(0)
  record_writer = RecordWriter(fastq_prefix)
  record_count = 0
  while True:
    record_1 = read_fastq_record(input_fastq)
    record_2 = read_fastq_record(input_fastq)
    record_count += 1
    #EOF will return None for a record
    if not record_1:
      if record_2:
        record_writer.write_unpaired_record(record_2)
      break
    elif not record_2:
      record_writer.write_unpaired_record(record_1)
      break
    if record_1.qname != record_2.qname:
      raise Exception("Input FASTQ not sorted. Paired records not adjacent")
    record_writer.write_paired_records(record_1, record_2)
    if not record_count % 10000000:
      print('Processed {} records'.format(str(record_count)))

def split_paired_fastq(fastq_1_path, fastq_2_path, fastq_prefix):
  """Converts a pair of FASTQ files to split chunks

  Args:
    fastq_1_path: Full path to first FASTQ file
    fastq_2_path: Full path to second FASTQ file
    fastq_prefix: Prefix for output FASTQ chunks
  """
  try:
    input_fastq_1 = gzip.open(fastq_1_path, 'r')
    input_fastq_2 = gzip.open(fastq_2_path, 'r')
    input_fastq_1.readline()
    input_fastq_2.readline()
  except IOError:
    input_fastq_1 = open(fastq_1_path, 'r')
    input_fastq_2 = open(fastq_2_path, 'r')
  input_fastq_1.seek(0)
  input_fastq_2.seek(0)
  record_writer = RecordWriter(fastq_prefix)
  record_count = 0
  while True:
    record_1 = read_fastq_record(input_fastq_1)
    record_2 = read_fastq_record(input_fastq_2)
    record_count += 1
    #EOF will return None for a record
    if not record_1:
      if record_2:
        record_writer.write_unpaired_record(record_2)
        continue
      break
    elif not record_2:
      record_writer.write_unpaired_record(record_1)
      continue
    if record_1.qname != record_2.qname:
      raise Exception("Input FASTQ not sorted. Paired records not adjacent")
    record_writer.write_paired_records(record_1, record_2)
    if not record_count % 10000000:
      print('Processed {} records'.format(str(record_count)))

def help_message():
  """Prints usage info"""
  print("\nUsage: python fastq_prep.py [out.fq prefix] [in.bam OR in.fq]")
  print("    OR python fastq_prep.py [out.fq prefix] [in_1.fq,in_2.fq]")
  print("\n*Multiple input files are separated by a comma")
  print("\nExample: python fastq_prep.py test_output test_input.bam")
  print("This command will generate test_output_XXXXX_R[1,2].fastq.gz\n")
  sys.exit()

def fastq_prep(output_prefix, input_files):
  """Main interface for converting files to split and compressed fastq files

  Determines type of input (Single BAM, Single FASTQ, or Paired FASTQ) and
  passes input files to correct conversion function.

  Args:
    output_prefix: full path prefix for output FASTQ chunks
    input_files: A list of input files to convert (1 or 2 files only)
  """
  print("\nConverting following file(s) to split and compressed FASTQ format:")
  for input_path in input_files:
    print("".join(["\t", input_path]))
  print("Output will be split into files with the following path structure:")
  print('\t{}_XXXXX_RX.fastq.gz'.format(output_prefix))
  if len(input_files) == 1:
    input_path = input_files[0]
    input_extension = os.path.splitext(input_path)[1].lower()
    if input_extension == ".bam" or input_extension == ".sam":
      split_bam(input_path, output_prefix)
    elif input_extension == ".fastq" or input_extension == ".fq":
      split_interleaved_fastq(input_path, output_prefix)
    else:
      raise Exception(" ".join(["Unknown file format:", input_extension]))
  elif len(input_files) == 2:
    input_path_1 = input_files[0]
    input_path_2 = input_files[1]
    split_paired_fastq(input_path_1, input_path_2, output_prefix)
  else:
    raise Exception('Incorrect number of input files given: {}'.\
      format(str(len(input_files))))
  print("Done")

if __name__ == "__main__":
  if len(sys.argv) != 3:
    help_message()
  fastq_prep(sys.argv[1], sys.argv[2].split(","))
