## Code base for FMR1 somatic instability analysis

### Usage:

#### Basic usage

```
$ python extract-reads.py -ref reference.fa -bed ./FMR1.bed -bam *.bam --analyse-methylation
```

#### Command line options
```
$ python extract_read.py -h

Extract the reads aligned at the FMR1 locus in a BAM file and analyse.

options:
  -h, --help            show this help message and exit
  -bam BAM [BAM ...]    Input BAM file from which the reads to be extracted
  -bed BED              Input regions file for which the reads are extracted
  -ref REF_FASTA        The reference fasta genome file
  -o OUTPUT, --output OUTPUT
                        Output file to write results. If not specified, prints to stdout.
  --aln-format ALN_FORMAT
                        Format of the alignment files. Choose from bam/cram/sam. Default: bam
  --plot                Plot the methylation levels of the reads. Default: False
  --analyse-methylation
                        Plot the methylation levels of the reads. Default: False
  --add-cigar           Plot the methylation levels of the reads. Default: False
  --haplotag HAPLOTAG   Haplotype tag to be used for the reads. Default: None
  --allele-bases        If the allele length should be reported in bases. By default allele length is reported in units. Default: False

```

### Citation


