from lib_ssw.pyssw import align_pair
from cigar_utils import convert_cigar, complement_cigar
import re


def reverse_complement(seq):
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]


def parse_modified_bases_from_tags(mm_tag, ml_tag, read):
    """
    Reconstructs pysam-like read.modified_bases dict
    from MM (modified bases) and ML (likelihoods) tags.

    Returns:
        { (base, mod_code, strand): [(position, likelihood), ...] }
    """

    if mm_tag is None:
        return {}

    # normalize ML to list
    if ml_tag is None:
        ml_vals = []
    elif isinstance(ml_tag, int):
        ml_vals = [ml_tag]
    else:
        ml_vals = list(ml_tag)

    ml_index = 0
    result = {}

    # MM has sections separated by ';'
    # each section encodes a modification type and delta positions
    # example: C+m?,0,5,3;
    sections = [s for s in mm_tag.split(';') if s]
    for section in sections:
        # parse header before commas: e.g., C+m?, C+h, A+a?
        if "," not in section: continue
        m = re.match(r'^([ACGTN])([+-])([A-Za-z?])(.*)$', section.replace(",", "", 1))
        if not m:
            continue
        base = section[0]                      # e.g., 'C'
        strand = section[1]                    # '+' or '-'
        mod_code = section[2]                  # e.g., 'm','h','?'

        strand = 0 if strand == '-' else 1
        key = (base, strand, mod_code)

        # extract all deltas after first comma
        # grab section starting after first comma
        deltas_part = section.split(',', 1)[1]
        deltas = [int(x) for x in deltas_part.split(',') if x != '']

        positions = []
        pos = 0  # deltas are relative; start -1 then add first delta

        ml_index = 0
        for d in deltas:
            likelihood = ml_vals[ml_index] if ml_index < len(ml_vals) else None
            ml_index += 1
            positions.append((pos, likelihood))

        result.setdefault(key, []).extend(positions)

    return result


def detect_flank(read, fasta, flank_len, repeat_start, repeat_end):
    """
        Identify the flanking sequences of the repeat region in the read
        Args:
            read: pysam read object
            fasta: pysam fasta object
            flank_len: length of the flanking sequence to be considered
            repeat_start: start coordinate of the repeat in reference
            repeat_end: end coordinate of the repeat in reference
        Returns:
            [start, end, sub_cigar]
            start and end coordinates of the flanking sequence in the read
            sub_cigar: CIGAR string of the repeat alignment to the reference in the read
    """
    print(f"Detecting flanks for read {read.query_name} at reference positions {repeat_start}-{repeat_end}")

    sub_cigar = ''

    # finding the upstream flank in the read 
    upstream = fasta.fetch(read.reference_name, repeat_start - flank_len, repeat_start)
    alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(read.query_sequence, upstream)

    # target_begin is the start of the alignment in the read
    # query_begin  is the start of the alignment in the reference flank
    # target_end   is the end of the alignment in the read
    # query_end    is the end of the alignment in the reference flank
    
    # if the alignment score is less than 90% of the flank length, return empty
    # making sure we identify the flanking sequence correctly in the read
    if alignment_score < int(2*(0.9*flank_len)): return [-1, -1, '']

    # reference position where the read (upstream flank) is starts to align
    read_reference_start = repeat_start - flank_len + query_begin - 1

    # the read alignment to the flank; the position where the read aligns to the start of the flank
    if target_begin > 1: sub_cigar += f'{target_begin - 1}S'

    # because we are aligning the flank (ref sequence) to the read
    # here we complement the cigar to get the alignment of the read to the reference flank sequence
    sub_cigar += complement_cigar(sCigar)

    # if the flank aligns completely this should be equal to the repeat start
    uflank_reference_end = read_reference_start + query_end - query_begin + 1
    
    # remaining bases in the flank that are not present in the read; these could be either deleted or substituted
    # these bases are at the end of the upstream flank
    uflank_remaining_length = flank_len - query_end
    uflank_read_end = target_end - 1
    # adding the remaining bases of the upstream flank as substitutions
    if uflank_remaining_length > 0:
        sub_cigar += f'{uflank_remaining_length}X'
        uflank_reference_end += uflank_remaining_length
        uflank_read_end += uflank_remaining_length

    # finding the downstream flank in the read
    downstream = fasta.fetch(read.reference_name, repeat_end, repeat_end + flank_len)
    alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(read.query_sequence, downstream)

    # target_begin is the start of the alignment in the read
    # query_begin  is the start of the alignment in the reference flank
    # target_end   is the end of the alignment in the read
    # query_end    is the end of the alignment in the reference flank

    # if the alignment score is less than 90% of the flank length, return empty
    # making sure we identify the flanking sequence correctly in the read
    if alignment_score < int(2*(0.9*flank_len)): return [-1, -1, '']
    
    dflank_cigar = ''
    
    # should be qual to the repeat end if the flank aligns completely
    dflank_reference_start = repeat_end + query_begin - 1
    
    # remaining bases in the flank that are not present in the read; these could be either deleted or substituted
    # these bases are at the start of downstream flank
    dflank_remaining_length = flank_len - query_end # number of bases in the flank that are not aligned
    dflank_read_start = target_begin - 1
    # adding the remaining bases of the downstream flank as substitutions
    if dflank_reference_start > repeat_end:
        dflank_cigar += f'{dflank_remaining_length}X'     # adding the remaining bases as substitutions
        dflank_reference_start -= (dflank_remaining_length)
        dflank_read_start -= (dflank_remaining_length)
    
    dflank_cigar += complement_cigar(sCigar)

    # instead of adding the remaining downstream and upstream flank bases as substitutions
    # they could be included with the repeat sequence to be aligned with the reference repeat sequence

    # read_repeat_sequence = read.query_sequence[uflank_read_end:dflank_read_start]
    start_idx = uflank_read_end + 1
    end_idx = dflank_read_start
    allele_length = end_idx - start_idx
    if (start_idx >= end_idx): return [-1, -1, '']

    # aligning the sequence between the two flanks to the reference repeat
    alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(read.query_sequence[start_idx: end_idx],
                                                                                                   fasta.fetch(read.reference_name, repeat_start, repeat_end))
    if target_begin - query_begin > 0: sub_cigar += f'{target_begin - query_begin}I'
    elif query_begin - target_begin > 0: sub_cigar += f'{query_begin - target_begin}D'
    if query_begin > 1: sub_cigar += f'{query_begin - 1}X'
    sub_cigar = complement_cigar(sCigar)
    remaining_subs = abs((repeat_end-repeat_start) - query_end)
    if remaining_subs > 0:
        sub_cigar += f'{remaining_subs}X'
        target_end += remaining_subs
        query_end += remaining_subs
    # if dflank_remaining_length > 0:
    #     sub_cigar += f'{dflank_remaining_length}D'
    if allele_length - target_end > 0:
        sub_cigar += f'{allele_length - target_end}I'

    # sub_cigar += dflank_cigar
    print(f"Read {read.query_name} - Upstream flank: {upstream}, Downstream flank: {downstream}, Sub CIGAR: {sub_cigar}")

    return [start_idx, end_idx, sub_cigar]


def parse_cigar(read, repeat_start, repeat_end, fasta, flank_len=50, motif_length=3):
    """
    Parses read alignment CIGAR and returns coordinates of the repeat within the read and the CIGAR
    corresponding to the repeat sequence.

    Args:
        cigar_tuples:   cigar tuples from pysam record
        read_start:     reference position of start of the read alignment
        repeat_start:   start coordinate of repeat in reference
        repeat_end:     end coordinate of repeat in reference
        ref_sequence:   reference sequence to which read is aligned
        query_sequence: sequence of the read

    Returns:
        start and end coordinates of read sequence aligning to the repeat region
        CIGAR string of alignment within repeat sequence
        [start, end, sub_cigar]
    """
    rpos = read.reference_start   # NOTE: The coordinates are 1 based in SAM; but pySAM uses 0 based
    qpos = 0            # starts from 0 the sub string the read sequence in python

    read_rpos = []
    # the read start and the repeat start are both on a 0 based coordinate system

    start_idx = -1; end_idx = -1
    sub_cigar = ''

    for c, cigar in enumerate(read.cigartuples):
        if cigar[0] == 4:   # soft clipped; do not affect the read position
           if rpos >= repeat_start - flank_len and rpos <= repeat_end + flank_len:
               # process the read to identify the flank
               start_idx, end_idx, sub_cigar = detect_flank(read, fasta, flank_len, repeat_start, repeat_end)
               return [start_idx, end_idx, sub_cigar]
           qpos += cigar[1]

        elif cigar[0] == 0 or cigar[0] == 7 or cigar[0] == 8: # match (both equals & difference)
            match_len = cigar[1]
            ctype = 'X' if cigar[0] == 8 else 'M'

            if start_idx == -1 and rpos + match_len > repeat_start:
                start_idx = qpos + (repeat_start - rpos)
                if rpos + match_len >= repeat_end:
                    end_idx = qpos + (repeat_end - rpos) 
                    sub_cigar += f'{repeat_end - repeat_start}{ctype}'
                else: sub_cigar += f'{rpos + match_len - repeat_start}{ctype}'

            elif start_idx != -1 and end_idx == -1:
                if rpos + match_len >= repeat_end:
                    end_idx = qpos + (repeat_end - rpos)
                    sub_cigar += f'{repeat_end - rpos}{ctype}'
                else:
                    sub_cigar += f'{match_len}{ctype}'

            # move both reference and query positions
            rpos += match_len; qpos += match_len

        elif cigar[0] == 1:     # insertion
            insert_length = cigar[1]

            if insert_length > motif_length and \
               (repeat_start - flank_len <= rpos <= repeat_start-1 or \
                repeat_end + 1 <= rpos <= repeat_end + flank_len):
                # process the read to identify the flanking sequences
                start_idx, end_idx, sub_cigar = detect_flank(read, fasta, flank_len, repeat_start, repeat_end)
                return [start_idx, end_idx, sub_cigar]

            if rpos == repeat_start and start_idx == -1:
                # if the insert is before the repeat include the sequence within the repeat
                start_idx = qpos
                sub_cigar += f'{insert_length}I'
            elif start_idx != -1 and end_idx == -1:
                sub_cigar += f'{insert_length}I'
            elif start_idx != -1 and rpos == repeat_end:
                sub_cigar += f'{insert_length}I'
                end_idx = qpos + insert_length
            # move the query position
            qpos += insert_length

        elif cigar[0] == 2:     # deletion
            # processing the deletion in the read is not necessary; need to check back
            deletion_length = cigar[1]
            if start_idx == -1 and rpos + deletion_length > repeat_start:
                start_idx = qpos
                if rpos + deletion_length >= repeat_end:
                    end_idx = qpos
                    sub_cigar += f'{repeat_end - repeat_start}D'
                else: sub_cigar += f'{rpos + deletion_length - repeat_start}D'

            elif start_idx != -1 and end_idx == -1:
                if rpos + deletion_length >= repeat_end:
                    end_idx = qpos
                    sub_cigar += f'{repeat_end-rpos}D'
                else: sub_cigar += f'{deletion_length}D'

            # move the reference position
            rpos += deletion_length

        if rpos > repeat_end:
            # if position moved beyond repeat
            return [start_idx, end_idx, sub_cigar]

    return [start_idx, end_idx, sub_cigar]


def get_perbase_methylation(read, start_idx, end_idx):
    """
    Get the per-base methylation level of the read.

    Args:
        read: pysam read object
        start_idx: start index of the repeat in the read
        end_idx: end index of the repeat in the read
    Returns:
        list of methylation levels for each base in the read
    """

    # Get methylation
    # dict (canonical base, strand, modification) -> list of tuples (position, quality)
    read_bms = [0]*(end_idx-start_idx); read_bmc = [0]*(end_idx-start_idx)
    called_bases = 0; ambiguous_bases = 0
    methylated_bases = 0; unmethylated_bases = 0

    if read.query_name == 'm84082_250131_064648_s1/116852188/ccs': # read exception with invalid MM tag
        return [read_bms, read_bmc, called_bases, ambiguous_bases, methylated_bases, unmethylated_bases]

    mods = read.modified_bases  # mods are a dict of base modifications
    if mods is None or len(mods) == 0:
        return [read_bms, read_bmc, called_bases, ambiguous_bases, methylated_bases, unmethylated_bases]

    if (('C', 0, 'm') not in mods) and (('C', 1, 'm') not in mods):
        return [read_bms, read_bmc, called_bases, ambiguous_bases, methylated_bases, unmethylated_bases]
    for modtype in mods:
        if modtype[0] == 'C' and modtype[2] == 'm':
            mods[modtype].sort(key=lambda x: x[0])
            for pos, qual in mods[modtype]:
                if pos >= start_idx and pos < end_idx: # Need to check the position logic
                    prob = (qual+1)/256
                    if prob > 0:
                        read_bms[pos - start_idx] = prob
                        if prob <= 0.33:
                            unmethylated_bases += 1; read_bmc[pos-start_idx] = 0.02
                        elif prob >= 0.66:
                            methylated_bases += 1; read_bmc[pos-start_idx] = 1
                        else:
                            ambiguous_bases += 1; read_bmc[pos-start_idx] = 0.5
                        called_bases += 1

    return [read_bms, read_bmc, called_bases, ambiguous_bases, methylated_bases, unmethylated_bases]
