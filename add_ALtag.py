import pysam
import sys

original_bam = sys.argv[1]
trgt_bam = sys.argv[2]
out_bam = sys.argv[3]

read_AL_tags = {}
trgt_bam = pysam.AlignmentFile(trgt_bam, "rb")
for read in trgt_bam.fetch(until_eof=True):
    if read.has_tag("AL"):
        read_AL_tags[read.query_name] = read.get_tag("AL")
trgt_bam.close()

original_bam = pysam.AlignmentFile(original_bam, "rb")
out_bam = pysam.AlignmentFile(out_bam, "wb", template=original_bam)
for read in original_bam.fetch(until_eof=True):
    if read.query_name in read_AL_tags:
        read.set_tag("AL", read_AL_tags[read.query_name])
    out_bam.write(read)
original_bam.close()
out_bam.close()