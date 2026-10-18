import pysam
import sys
import warnings

# At the moment, is currently redundant as Vardict Automatically handles this
# However, will be re-incorporated for the SV breakpoint detection, to be implemented in the future
# Even if a read poorly maps, if its mate or supplementary alignment align well, it will be kept
# Ensuring that the breakpoints could be identified

def keep_read(read):
    mismatches = 0
    # intronic_indel = False

    # list of tuples, consisting of 0 based offset from the start of the read sequence followed by the 0-base reference position
    # with seq returns the additional reference sequence,

    aligned_pairs = read.get_aligned_pairs(with_seq=True)
    left_clip = 0
    right_clip = 0

    # note that read.get_aligned_pairs(with_seq=True) also gets soft clipped read bases, however they get displayed like this:
    # for example, for the very first read base
    # (0, None, None)
    # so 7S144M appears as (0, None, None), (1, None, None), (2, None, None), (3, None, None), (4, None, None), (5, None, None), (6, None, None)

    read_length = read.query_length

    # an op value of "0" means "M" in the cigar
    # an op value of "4" means "S" in the cigar
    if read.cigartuples:
        op, length = read.cigartuples[0]
        if op == 4:
            left_clip = length
        
        op, length = read.cigartuples[-1]
        if op == 4:
            right_clip = length

    for read_i, ref_pos, ref_base in aligned_pairs:
        # For inserts, deletions, skipping either query or reference position may be None.
        # For padding in the reference, the reference position will always be None.
        # do not assess indels, just assess based on direct base mismatches
        # large indels could easily build up exon/intron mismatches, when the event could be real
        if read_i is not None and ref_base is not None:
            # identifying if current base is in a soft clipped region
            if read_i < left_clip or read_i >= read_length - right_clip:
                continue
            read_base = read.query_sequence[read_i]
            if ref_base.upper() != read_base.upper():
                mismatches += 1

    # note that increasing the exon_mismatches + intron_mismatches (while disablnig other checs) 
    # will not save many reads
    # for sample "13"
    # current value results in 206 singletons
    # increasing it to 10 results in 150 singletons
    # increasing it to at least 20 results in 136 singletons
    # with all the below checks, resuts in 218 singletons

    if mismatches > 5:
        return False

    return True

input_bam = sys.argv[1]
output_bam = sys.argv[2]

in_bam_reader = pysam.AlignmentFile(input_bam, "rb")

out_bam_writer = pysam.AlignmentFile(output_bam, "wb", template=in_bam_reader)

num_secondary_or_supplementary = 0

num_processed_records = 0
for read in in_bam_reader.fetch("chrX"):
    num_processed_records += 1
    if read.has_tag("MD") is False:
        assert read.is_unmapped is True # unmapped read usually given the coordinates of its mapped mate
        continue
    try:
        if keep_read(read):
            out_bam_writer.write(read)
    except Exception as e:
        warnings.warn(e)
        warnings.warn("An error arose when assessing: ${s}, read will be skipped", read.query_name)
    if read.is_secondary or read.is_supplementary:
        num_secondary_or_supplementary += 1

print("Processed %{s} records", num_processed_records)
print("There were " + str(num_secondary_or_supplementary) + " reads that were either secondary or supplementary.")

in_bam_reader.close()
out_bam_writer.close()
