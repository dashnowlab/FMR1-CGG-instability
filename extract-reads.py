import pysam
import os, sys, argparse
import statistics as stats
import numpy as np

from data_viz import *
from process_read import *
# from haplotype_reads import cluster_reads


"""
NOTE: Though samtools uses a position with a 1 based coordinate system. read.reference_start returns the position w.r.t 
      a 0-based coordinate system.
"""

def parse_args():
    parser = argparse.ArgumentParser(prog='extract-reads', description="Extract the reads aligned at the FMR1 locus in a BAM file and analyse.")

    parser.add_argument('-bam', required=True, nargs="+", type=str, dest="bam", help="Input BAM file from which the reads to be extracted")
    parser.add_argument('-bed', required=True, type=str, dest="bed", help="Input regions file for which the reads are extracted")
    parser.add_argument('-ref',  required=True, type=str, dest="ref_fasta", help="The reference fasta genome file")

    parser.add_argument('-o', '--output', type=str, default=None, dest="output", help="Output file to write results. If not specified, prints to stdout.")

    parser.add_argument('--aln-format', default='bam', type=str, dest="aln_format", help="Format of the alignment files. Choose from bam/cram/sam. Default: bam")
    parser.add_argument('--plot', action='store_true', dest="plot", help="Plot the methylation levels of the reads. Default: False")
    parser.add_argument('--analyse-methylation', action='store_true', dest="analyse_methylation", help="Plot the methylation levels of the reads. Default: False")
    parser.add_argument('--add-cigar', action='store_true', dest="add_cigar", help="Plot the methylation levels of the reads. Default: False")

    parser.add_argument('--haplotag', default=None, type=str, dest="haplotag", help="Haplotype tag to be used for the reads. Default: None")
    parser.add_argument('--allele-bases', action='store_true', dest="allele_bases", help="If the allele length should be reported in bases. By default allele length is reported in units. Default: False")
    args = parser.parse_args()

    return args


def count_mm_positions(mm_string):
    # crude but works: count commas
    return mm_string.count(",") if mm_string else 0


def bam_check_tags(bam_file, format, args):
    """
    Check if the BAM files have the required tags for ATaRVa.
    
    Args:
        bam_file (str): Path to the BAM file.
        parameters (argparse.Namespace): Parsed command line arguments.
    
    Returns:
        None
    Raises:
        ValueError: If the BAM file does not contain the required tags.
    """
    print(f"Checking alignment tags in {os.path.basename(bam_file)}", file=sys.stderr)
    sample_size = 1000 # sample size of reads to check for tags
    
    reads_sampled = 0 # sample set of reads to look for tags
    if format == 'rc':
        aln_file = pysam.AlignmentFile(bam_file, format, reference_filename=args.ref_fasta, check_sq=False)
    else:
        aln_file = pysam.AlignmentFile(bam_file, format, check_sq=False)
    cs_tag = False; md_tag = False; cigar_tag = False
    
    for read in aln_file.fetch():
        # 0x400 - read is PCR or optical duplicate
        # 0x100 - not primary alignment
        if (read.flag & 0X400) or (read.flag & 0X100): continue 
        reads_sampled += 1
        cigar = read.cigarstring
        
        if read.has_tag('cs') or read.has_tag('CS'):
            cs_tag = True
        
        elif (cigar!=None) and (('X' in cigar) or ('=' in cigar)):
            cigar_tag = True
        
        elif read.has_tag('MD'):
            md_tag = True
        
        if reads_sampled >= sample_size:
            break

    aln_file.close()

    return [cs_tag, md_tag, cigar_tag]


def extract_reads(args, aln_format):
    """
    Given the set of repeat loci. The function extracts reads aligning to the repeat locus
    within a given set of BAMs

    Args:
        bed_file:   list of input of repeat regions as a BED file
        bam_files:  list of input BAM files (Ribbit input)
        aln_format: alignment file format (bam/sam/cram)
    """

    # stores the alignment tags present for each BAM file
    # key: bam_id, value: [cs_tag, md_tag, cigar_tag]
    bam_aln_tags = {}
    for bam_file in args.bam:
        try:
            bam_id = os.path.basename(bam_file).decode('utf-8').split('.')[0]
        except AttributeError:
            bam_id = os.path.basename(bam_file).split('.')[0]
        if bam_id not in bam_aln_tags:
            bam_aln_tags[bam_id] = bam_check_tags(bam_file, aln_format, args)

    # bams = [pysam.AlignmentFile(bam_file, aln_format, check_sq=False) for bam_file in args.bam]
    bams = []
    for bam_file in args.bam:
        if aln_format == 'rc':
            bam = pysam.AlignmentFile(bam_file, aln_format, reference_filename=args.ref_fasta, check_sq=False)
        else:
            bam = pysam.AlignmentFile(bam_file, aln_format, check_sq=False)
        bams.append(bam)
    fasta = pysam.FastaFile(args.ref_fasta)

    if args.output is not None:
        sys.stdout = open(args.output, 'w')

    flank_len = 50

    read_qualities = []
    read_methylation = []

    bam_names = {}
    with open('./sample-metadata.tsv') as fh:
        header = fh.readline().strip().split('\t')
        for line in fh:
            line = line.strip().split('\t')
            bam_names[line[-2]] = f'{line[2]}'

    header = ['#bam_id', 'loc_id', 'read_id', 'read_repeat_start', 'read_repeat_end', 'allele_length', 'base_qual', 'trgt_haplotype', 'direction']
    if args.add_cigar: header += ['cigar']
    if args.haplotag is not None: header += ['haplotype']
    if args.analyse_methylation:
          header += ['median_meth', 'called_bases', 'methylated_bases', 'umethylated_bases', 'ambiguous_bases']
    print(*header, sep='\t')

    with open(args.bed) as fh:
        for line in fh:
            line = line.strip().split('\t')
            chrom = line[0]
            repeat_start = int(line[1]); repeat_end = int(line[2])
            motif_len = len(line[3])

            sample_allele_counts = {}

            for bam in bams:
                bam_file = bam.filename
                # Getting the base file name
                try:
                    bam_id = os.path.basename(bam_file).decode('utf-8').split('.')[0]
                except AttributeError:
                    bam_id = os.path.basename(bam_file).split('.')[0]
                if bam_id not in sample_allele_counts:
                    sample_allele_counts[bam_id] = {}
                if chrom not in bam.references: continue

                # clustered_reads = cluster_reads(bam, chrom, repeat_start-flank_len, repeat_end+flank_len, bam_aln_tags[bam_id])

                reads = bam.fetch(chrom, repeat_start-flank_len, repeat_end+flank_len)
                read_data = []; sample_reads_bqs = []; sample_reads_bms = []

                for read in reads:
                    if read.is_unmapped or read.is_secondary or read.is_supplementary: continue

                    start_idx = -1; end_idx = -1

                    # adjust the potential start and end of the read based on the soft-clipped bases
                    reference_start = read.reference_start - read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else read.reference_start
                    reference_end   = read.reference_end + read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else read.reference_end

                    # read.reference_start is 0 based; though the position in the SAM file is 1 based
                    if reference_start < repeat_start-flank_len and reference_end > repeat_end + flank_len:
                        start_idx, end_idx, sub_cigar = parse_cigar(read, repeat_start, repeat_end, fasta, flank_len, motif_len)

                    # if start_idx and end_idx are not found in any reads
                    if start_idx == -1 or end_idx == -1:
                        continue

                    # Get base call qualities
                    base_qualities = read.query_qualities[start_idx:end_idx]
                    if len(base_qualities) == 0: continue

                    if np.mean(base_qualities) < 25: continue  # filter out low quality reads

                    if args.allele_bases:
                        allele_len = (end_idx - start_idx)
                    else:
                        allele_len = round((end_idx - start_idx)/motif_len, 2)
                    # allele_len = (end_idx - start_idx)

                    read_data.append([bam_id, f"{chrom}:{repeat_start}-{repeat_end}", read.query_name, start_idx, end_idx,
                                      allele_len, round(np.mean(read.query_qualities[start_idx:end_idx]), 2)])
                    if args.add_cigar: read_data[-1].append(sub_cigar)

                    # Get AL tag added by TRGT
                    if read.has_tag('AL'):
                        read_data[-1].append(read.get_tag('AL'))
                    else: read_data[-1].append('NA')

                    if read.is_reverse: read_data[-1].append('-')
                    else: read_data[-1].append('+')

                    read_data[-1].append(read.query_sequence[start_idx:end_idx])
                    sample_reads_bqs.append((list(base_qualities), read_data[-1][-1]))

                    if args.analyse_methylation and (read.has_tag('MM') and read.has_tag('ML')) and len(read.get_tag('ML')) > 0:
                        # Get methylation
                        read_bms, read_bmc, called_bases, ambiguous_bases, methylated_bases, unmethylated_bases = get_perbase_methylation(read, start_idx, end_idx)

                        sample_reads_bms.append(read_bmc)
                        # median meth
                        read_bms = [meth for meth, base in zip(read_bms, read.query_sequence[start_idx:end_idx]) if base == 'C']
                        if len(read_bms) > 0: med_meth = round(stats.median(read_bms), 2)
                        else: med_meth = None
                        read_data[-1] += [med_meth]
                        read_data[-1].append(called_bases)
                        if called_bases > 0:
                            read_data[-1].append(methylated_bases)
                            read_data[-1].append(unmethylated_bases)
                            read_data[-1].append(ambiguous_bases)
                        else:
                            read_data[-1].append(None)
                            read_data[-1].append(None)
                            read_data[-1].append(None)

                    # adding the read sequence to the read_data
                    # read_data[-1].append(read.query_sequence[start_idx:end_idx])
                    if args.haplotag is not None:
                        tags = [x[0] for x in read.tags]
                        if args.haplotag in tags:
                            read_data[-1].append(read.get_tag(args.haplotag))
                        else:
                            read_data[-1].append('NA')
                for rd in read_data: print(*rd, sep='\t', file=sys.stdout)

                if args.plot:
                    indices = [i for i, v in sorted(enumerate(sample_reads_bqs), key=lambda x: len(x[1][0]))]
                    indices = sorted(range(len(sample_reads_bqs)), key=lambda i: len(sample_reads_bqs[i][0]))
                    sample_reads_bqs = [sample_reads_bqs[i][0] for i in indices]
                    sample_reads_bms = [sample_reads_bms[i] for i in indices]
                    plot_bqsvsbms(sample_reads_bqs, sample_reads_bms, bam_id, f'../plots/{"-".join(bam_names[bam_id].split(" "))}')

    for bam in bams: bam.close()
    fasta.close()
    if args.output is not None: sys.stdout.close()


if __name__ == "__main__":
    args = parse_args()

    aln_format = 'rb'
    if args.aln_format == 'cram': aln_format = 'rc'
    elif args.aln_format == 'sam': aln_format = 'r'

    extract_reads(args, aln_format)
