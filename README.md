# RNA-seq Unannotated Transcripts Detection Pipeline

Scripts for detecting unannotated transcripts from an RNA-seq library using a comprehensive pipeline. The pipeline is using `find_unannotated_transcripts.sh`, which uses various other scripts for different stages of the process. Below is an overview of the pipeline and instructions for usage.

## Pipeline Overview

### Step 1: Read Extraction
1. **Extract Genome Coverage Reads**: bedtools genomecov extracts the number of reads for every base in the genome.

### Step 2: Initial Transcript Identification
2. **Identify Transcripts with Minimum Coverage**: 
   - **2a**: Utilizing `genomecov_transcripts.py`, regions with a minimum specified number of reads per base (e.g., 100 reads) are extracted, allowing reads to drop for a specified number of bases (e.g., 10 bases).
   - **2b**: Results shorter than 40 base pairs (bps) are removed.

### Step 3: Reverse Transcript Identification
3. **Reverse Transcript Identification**:
   - **3a**: Repeat Step 2a in reverse order, from the last base to the first base.
   - **3b**: Remove results shorter than 40 bps.

### Step 4: Transcript Merging
4. **Merge Transcripts**: Using bedtools merge, the transcripts resulted from Steps 2 and 3 are merged.

### Step 5: Annotation Filtering
5. **Filter Annotated Genes**:
   - **5a**: Known annotated genes are removed from the results, leaving only unannotated regions.
   - **5b**: Results shorter than 40 bps are filtered out.

### Step 6: Refinement
6. **Refinement**: Repeat Steps 1-4 on the unannotated regions as the "genome" to improve accuracy.

### Step 7: Sequence Complexity Filtering
7. **Filter Low Complexity Sequences**: Utilize `prinseq-lite` to filter out sequences with low complexity (e.g., "AAAAA...", "AGAGAG...", "GTAGTAGTA...").

## Usage
1. Clone the repository to your local machine.
2. Ensure all dependencies are installed:
   - bedtools
   - prinseq-lite
3. Run `find_unannotated_transcripts.sh` <bam file> <known genes annotations bed file> [max reads (default: 100)] [max gaps (default: 10)]
4. Further analysis can be performed on the resulting unannotated transcripts using tools like CPC2, LncFinder, etc.

## Requirements
- Python (for `extract_genomecov_reads.py` and `genomecov_transcripts.py`)
- bedtools
- prinseq-lite

## Note
- Adjust parameters in the scripts according to your specific requirements, such as read coverage thresholds and sequence length filters.
