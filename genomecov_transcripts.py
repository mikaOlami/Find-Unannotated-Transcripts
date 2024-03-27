######################################
# Calculate Transcipts Using Genomecov
# Author: Mika Olami
# Date: 25/04/2021
######################################
# !/usr/bin/env python
import sys
import math

# Initialize variables
current_start = 0
current_end = 0
current_chr = ''
cov_peak = 0
bed_dict = {}
gap_count = 0

# Check command line arguments
if len(sys.argv) < 4:
    exit("\nUsage: python3 get_genomecov_transcripts.py [genomecov results] "
         "[min coverage # per transcript] [max gap length]\n")

# Parse command line arguments
genomecov_file = sys.argv[1]
min_coverage = int(sys.argv[2])
max_gap = int(sys.argv[3])

# Read the genomecov results
with open(genomecov_file, "r") as genomecov_coverage:
    # Iterate through each line
    for base in genomecov_coverage:
        chrom, bp, cov = base.strip().split("\t")
        bp = int(bp)
        cov = float(cov)

        # Check if the bases are sequential
        if current_end + 1 != bp and current_end - 1 != bp:
            # If there is a jump, split transcripts and reset counters
            if len(current_chr) > 0:
                if current_chr not in bed_dict:
                    bed_dict[current_chr] = []
                if current_start != current_end:
                    bed_dict[current_chr].append([current_start, current_end])
                current_chr = ''
                gap_count = 0
                cov_peak = 0
                min_coverage = int(sys.argv[2])
                current_start = 0
                current_end = 0

        # Check coverage
        if cov >= min_coverage:
            cov_peak = max(cov, cov_peak)
            min_coverage = math.ceil(cov_peak * 0.05)

            # Check chromosome
            if len(current_chr) == 0 or current_chr != chrom:
                if len(current_chr) > 0 and [current_start, current_end] not in bed_dict[current_chr]:
                    bed_dict[current_chr].append([current_start, current_end])
                current_chr = chrom
                cov_peak = cov
                current_start = bp

            current_end = bp
            gap_count = 0  # Reset gap counter
            continue
        elif len(current_chr) > 0 and gap_count < max_gap:
            # Allowed gaps
            gap_count += 1
            continue
        else:
            # Passed the max gap length
            gap_count += 1
            if gap_count > max_gap:
                if len(current_chr) > 0:
                    if current_chr not in bed_dict:
                        bed_dict[current_chr] = []
                    if current_start != current_end:
                        bed_dict[current_chr].append([current_start, current_end])
                    current_chr = ''
                    gap_count = 0
                    cov_peak = 0
                    min_coverage = int(sys.argv[2])

# Check for any remaining transcripts
if len(current_chr) > 0:
    if current_chr not in bed_dict:
        bed_dict[current_chr] = []
    if current_start != current_end:
        bed_dict[current_chr].append([current_start, current_end])

# Write the output to a BED file (chr   start   end)
for key in sorted(bed_dict.keys()):
    for line in sorted(bed_dict[key]):
        print("%s\t%d\t%d" % (key, line[0], line[1]))
