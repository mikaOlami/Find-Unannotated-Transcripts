#########################################################################################
# The complete bash pipeline to extract unannotated transcripts from a sequence library #
#########################################################################################

# Function to display usage information
usage() {
    echo;
    echo "Usage: $0 <sorted_bam> <known_annots> [min_coverage] [max_gap]";
    echo "Arguments:";
    echo "  <sorted_bam>: Sorted library BAM file";
    echo "  <known_annots>: Genome's known annotations (BED file)";
    echo "  [min_coverage]: Minimum reads for transcript detection (Default: 100)";
    echo "  [max_gap]: Maximum gaps allowed for transcript detection (Default: 10)";
    echo;
    exit 1;
}

# Check if the required number of arguments are provided
if [ "$#" -lt 2 ]; then
    usage;
fi


sorted_bam=$1;    # script input - sorted library bam file
known_annots=$2;  # genome's known annotations (bed file)
base_file=$(basename "$sorted_bam" | sed 's/.bam//g');  # library base file name
genomecov="$base_file".genomecov;  # genomecov file name
genomecov_rev="$genomecov".rev;    # reverse genomecov file name
min_coverage="${3:-100}";          # minimum reads for transcript detection (Default: 100)
max_gap="${4:-10}";                # max gaps allowed for transcript detection (Default: 10)


transcripts_front="$base_file".front.bed;
transcripts_rev="$base_file".rev.bed;
transcripts_merged="$base_file".merged.bed;
no_overlaps="$base_file".no_overlaps.bed;

genomecov_no_overlaps="$genomecov".noOverlaps;
genomecov_no_overlaps_rev="$genomecov_no_overlaps".rev;

transcripts_filt_front="$base_file".filtered.front.bed;
transcripts_filt_rev="$base_file".filtered.rev.bed;
transcripts_filt_merged="$base_file".filtered.bed;

pyscript_path=".";  # edit the location of the python scripts

#############################################################

# Running genomecov
bedtools genomecov -ibam "$sorted_bam" -d > "$genomecov";

# Reversing genomecov file
tac "$genomecov" > "$genomecov_rev";

# Calculate genomecov transcripts, and filter results shorter than 40bps
python3 "$pyscript_path"/genomecov_transcripts.py "$genomecov" "$min_coverage" "$max_gap" | awk '$3-$2>=40' > "$transcripts_front";

# Calculate reverse genomecov transcripts, and filter results shorter than 40bps
python3 "$pyscript_path"/genomecov_transcripts.py  "$genomecov_rev" "$min_coverage" "$max_gap" | awk '{print $1 "\t" $3 "\t" $2}' | awk '$3-$2>=40' > "$transcripts_rev" ;

# Merge front and reverse transcripts
bedtools merge -i <(cat "$transcripts_front" "$transcripts_rev"  | sort -k1,1 -k2,2n) > "$transcripts_merged";

# Remove the annotated regions from the results, and filter transcripts shorter than 40bps
bedtools subtract -a "$transcripts_merged" -b "$known_annots" | awk '$3-$2>=40' > "$no_overlaps";

# Extract the reads of the non-overlapping transcripts
python3 "$pyscript_path"/extract_genomecov_reads.py "$genomecov" "$no_overlaps" > "$genomecov_no_overlaps";

# Reversing the genomecov of the non-overlapping transcripts
tac "$genomecov_no_overlaps" > "$genomecov_no_overlaps_rev";

# Second round of the python script on the previous results
python3 "$pyscript_path"/genomecov_transcripts.py "$genomecov_no_overlaps" "$min_coverage" "$max_gap" | awk '$3-$2>=40' > "$transcripts_filt_front" ;

# Second round of the python script on the previous results (reversed)
python3 "$pyscript_path"/genomecov_transcripts.py "$genomecov_no_overlaps_rev"  "$min_coverage" "$max_gap" | awk '{print $1 "\t" $3 "\t" $2}' | awk '$3-$2>=40' > "$transcripts_filt_rev";

# Merge filtered front and reverse transcripts
bedtools merge -i <(cat "$transcripts_filt_front" "$transcripts_filt_rev"  | sort -k1,1 -k2,2n) > "$transcripts_filt_merged";

# Delete temporary files
rm "$transcripts_front" "$transcripts_rev" "$transcripts_merged" "$no_overlaps" "$genomecov_no_overlaps" "$genomecov_no_overlaps_rev" "$transcripts_filt_front" "$transcripts_filt_rev";

# Done