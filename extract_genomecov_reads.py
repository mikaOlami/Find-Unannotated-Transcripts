#!/usr/bin/env python

import sys

def extract_coverage(genomecov_file, bed_file):
    # Create a dictionary to store coverage values for each position
    coverage_dict = {}

    # Read the genomecov file
    with open(genomecov_file, "r") as genomecov_coverage:
        for line in genomecov_coverage:
            chrom, pos, cov = line.strip().split("\t")
            pos = int(pos)
            cov = int(cov)
            coverage_dict.setdefault(chrom, {})[pos] = cov

    # Read the bed file and print the coverage for each nucleotide
    with open(bed_file, "r") as bed_coordinates:
        for line in bed_coordinates:
            line = line.strip().split("\t")
            chrom = line[0]
            start = int(line[1])
            end = int(line[2]) 

            # Print coverage for each nucleotide in the specified range
            for pos in range(start, end + 1):
                if chrom in coverage_dict and pos in coverage_dict[chrom]:
                    print(f"{chrom}\t{pos}\t{coverage_dict[chrom][pos]}")
                else:
                    # Print 0 coverage if the position is not found in the genomecov file
                    print(f"{chrom}\t{pos}\t0")

if __name__ == "__main__":
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python3 extract_coverage.py <genomecov_file> <bed_file>")
        sys.exit(1)

    genomecov_file = sys.argv[1]
    bed_file = sys.argv[2]

    # Call the function to extract coverage
    extract_coverage(genomecov_file, bed_file)
