# Meta-MARC
## Metagenomic Markov models for Antimicrobial Resistance Characterization
Created by the Microbial Ecology Group at Colorado State University

Developer: Steven Lakin (Steven.Lakin@colostate.edu)

## Description
Meta-MARC is a set of profile Hidden Markov Models develoepd for the purpose of screening and profiling
resistance genes in DNA-based metagenomic data.  This tool was developed for the characterization of various
resistance classes, mechanisms, and gene/operon groups from raw sequencing data much in the way that
microbiome tools profile the bacterial taxonomy of metagenomic samples.  Meta-MARC is not intended to be
used as a final annotation or all-in-one tool; this software simply offers the user another view of complex 
metagenomic data.  Users are encouraged to perform further statistical testing
on their data to validate results obtained from Meta-MARC.

## Dependencies
The following dependencies are for use of the mmarc workflow scripts.  If you wish to use your own workflow,
then copy the HMM files in meta-marc/src/HMMs to the desired location on your server/computer for
use with your installation of HMMER.

**HMMER**

**Python 3.4** or higher
 - matplotlib Python package
 - numpy Python package

**Bash shell and the GNU coreutils**
 - workflow scripts were tested on GNU Bash v4.3.11 and later
 - scripts were tested on Ubuntu and Red Hat Linux (no testing was performed for other operating systems)

## Installation (Recommended)
The Meta-MARC source code available here contains the pre-compiled Hideen Markov Models (HMMs) for use with
the HMMER suite of tools.  The **Build.sh** script in the meta-marc directory will automatically download
the appropraite HMMER binaries and create the necessary symlinks for Meta-MARC to run.
If you would like to use your own HMMER installation, simply skip to the Usage section and the mmarc
executable will detect HMMER in your PATH.

Note: If you use the Meta-MARC HMM files independently of the mmarc workflow, you will need to annotate your
own HMMER output.  **The annotation file for the Meta-MARC models is located here: meta-marc/src/mmarc\_model\_annotations.tsv**

Steps for installation:

Either download the zip file from GitHub and extract, or use the following command to clone into the repository.  
Then cd into the directory and run Build.sh.  Answer the prompts as necessary.
```
git clone --single-branch https://github.com/lakinsm/meta-marc.git
cd meta-marc
./Build.sh
```

## Meta-MARC structure
The primary bash workflow script **mmarc** is located in the bin file and uses relative pathing.  Either run it
directly from the meta-marc/bin file or use a symlink.

**Recommended options**

Meta-MARC was developed for running on FASTQ data or assembled contigs of metagenomic origin.  The pipeline does not account for base
 quality scores or orientation of reads.  For unassembled reads in FASTA or FASTQ format, it is recommended to use the 
deduplication (-d) option.  The user can decide how reads are counted: the multihit correction
option keeps a 1:1 ratio of reads in the input file to hits on the output file (recommended).  Any read that hits multiple
models below the desired E-value threshold will be split amongst the models hit when the -m option is enabled.
When the -m option is not enabled, any read hitting to multiple models below the E-value threshold will be
counted multiple times.  Best results were achieved in the publication using the multihit correction option.

For assembled reads in FASTA format, it is recommended to not use multihit correction, as the pipeline
does not take into account how many reads were used to assemble a given contig region.  Additionally,
deduplication of the input file will not likely have an effect on the output with assembled data.

See the publication for a description of results obtained using the FASTQ and contig pipelines.

**K-mer screen description**

WARNING: Not validated and not included in the publication

The k-mer screen hashes forward and reverse k-mers from the input database and runs a sliding window along
each read from the input analysis FASTA/FASTQ file.  If a single k-mer matches (exactly) to the hashtable
of reference k-mers, the read is output and run through the HMMER pipeline.  To see the extent to which
the k-mer screen affects sensitivity and specificity of the reuslting calls by HMMER, see the publication.

The k-mer screen saves substantial time, especially on large input files where many reads do not match to
the HMMs.  As optimized as HMMER is, its pipeline is still computationally costly, and the k-mer screen
helps to reduce the overall runtime significantly.

**Model set description**

The Meta-MARC models are divided into three groups: groupI, groupII, and groupIII.

*groupI: 262 models*

GroupI contains a core set of models and was constructed using only sequences from curated databases with known or annotated
resistance mechanisms.  This includes the Comprehensive Antibiotic Resistance Database (CARD) and the
database of beta-lacatamase sequences maintained by the Lahey clinic.  Additionally, only sequences that
clustered together from a sequence distance cutoff of 80% using USEARCH were included.  Mechanisms
requiring the presence of a specific SNP were not included in groupI (such as elfamycin EF-tu and
trimethoprim dihydrofolate reductase resistance).

*groupII: 316 models*

GroupII includes all of groupI and models containing mechanisms based on one or more SNPs.
Running this set of models with the Meta-MARC workflow script will automatically screen the HMMER output
for the SNPs described in the literature for HMMER hits overlapping that location in the groupII models.
Reads without described resistance SNPs will not be counted as resistance genes and will be output as a
separate file.

*groupIII: 678 models*

GroupIII includes all of groupII and the addition of sequences from groupsI-II that did not cluster
with any other sequences (singletons).  These singleton genes were augmented with other sequences using an
80% sequence similarity threshold with BLASTn against the NCBI nt database.  Sequences were reclustered with
the augmented sequences and these resulting clusters were used to construct the additional models for groupIII.
All models that are augmented with BLAST sequences in groupIII are annotated as such with "blast_augmented"
at the end of the model name.

**Annotation description**

The Meta-MARC models have been annotated in a hierarchical, four-level manner:

 - Class: annotation at the class level describes the class of antibiotics to which this entry confers resistance (e.g. tetracyclines, betalactams)
 - Mechanism: annotation at the mechanism level describes the molecular function that confers resistance to a given antibiotic class (e.g. Class D betalactamases)
 - Group: annotation at the group level describes the gene or operon group under which this entry falls (e.g. Aph-6'', OXA betalactamase)
 - HMM: this level describes the name of the Meta-MARC model, constructed from sequences of at least 80% similarity
 
Annotations for all model sets can be found in the master annotation file in meta-marc/src/mmarc\_model\_annotations.tsv

## Meta-MARC workflow

Example for FASTQ (compressed or uncompressed) or FASTA files of reads:
```
mmarc -i input_file.fastq.gz -o output_directory -f output_filename_no_extension -d -l 3 -m -t 4
```

Example for paired end FASTQ files (compressed or uncompressed):
```
mmarc --i1 foward_file.fastq.gz --i2 reverse_file.fastq.gz -o output_directory -f output_filename -d -l 3 -m -t 4
```

Example for FASTA files containing assembled contigs:
```
mmarc -i input_contig_file.fasta -o output_directory -f output_filename_no_extension -l 3 -t 4
```

**Detailed description of options**
```
-d | --dedup
De-duplicate the fastq before running HMMER; this saves significant time, since HMMER's screening algorithms
are relatively costly compared to detecting duplicate reads in a FASTQ file.  The de-duplicated reads will
automatically be re-inflated (the counts will be returned to normal) by the mmparse.py script after HMMER
completes its run.  It is recommended to use this in all cases where raw reads are being analyzed.

-e | --evalue	FLOAT
The E-value under which results from HMMER will be kept in scientific form (1e-3); all HMMER hits above this value will be ignored.
The default setting is 10, which essentially includes most/all HMMER output.  Change this if you would like
a more stringent threshold for detection.  Note that in our publication describing Meta-MARC, the best
performance for general profiling of the resistome was achieved without filtering on E-value.

-f | --filename	STR
The output base filename to be used for the intermediate file stages in the mmarc workflow (no file extension)

-g | --graphs	DIR
Directory to output graphs of HMM coverage (a coverage graph); ONLY USE WITH -w ENABLED.  Warning: there will
be potentially as many output graphs are there are models in the model set (at most 678).

-i | --input	FILE
The input file in FASTA or FASTQ format.  See the above recommended options for more information.

   --i1 FILE
   --i2 FILE
   Use the above options for paired-end input (this will only work for fastq or gzip'd fastq files).

-k | --kmer	INT**
The length of k-mer used in the k-mer screen.  If you choose to use this option (not recommended),
it is recommended is to use 16- or 15-mers based on sensitivity analysis.

-l | --level	INT
The number describing the model set to use; choices are 1, 2, or 3.  The models are described in detail
in the publication and above in the "Model set description" section.

-m | --multicorrect
This is a flag for whether or not to correct for multiple hits per read.  For example, if a read from a
FASTQ file hits to more than one model in the model set, the hit (1 hit) will be divided across those models
(1 read hitting to 3 models results in 0.333 counts for each model).  This option retains a 1:1 mapping
from reads to output counts.  See the recommended options section for suggestions on how to use this option.

-o | --output	DIR
The top-level directory for output to be written.  Each sample run through the workflow script will generate
1-4 files depending on the options used.  The basename described in the -f option will be used for the base
filename at each step.

-t | --threads	INT
Number of threads to use for the k-mer screen and HMMER pipeline.  The default is 1.

-w | --skewness	FILE
Run a "skewness" analysis and output the results to this file.  The -g option must be used when -w is
enabled; graphs of coverage for each model will be generated.  Metrics calculated include a percentage of
maximum Shannon's entropy and a percentage of closeness to the n-sphere uniform distribution using the
2-norm.  This option is useful for seeing how reads are distributed across each model.
```


## Rebuilding the Meta-MARC HMMs

For users who wish to rebuild the Meta-MARC HMM files from the original FASTAs, this can be done using the
mmarc\_build\_models.sh script under the meta-marc/src/building directory on the building branch.  
**It is the user's responsibility to have the necessary dependencies installed to rebuild the HMM files.**
This includes obtaining a binary file for USEARCH (the models were originally built with USEARCH v8.1.1861_i86linux32).
USEARCH requires an end-user license agreement, and therefore is not included by default in the Meta-MARC files.
The binary file must be renamed to "usearch" and be placed in the same directory as the building script.


