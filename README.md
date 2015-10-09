# poremap version 0.1.1
Pipeline for nanopore read alignment and statistics
```
Usage:

  poremap.py [-h] --bindir DIR [--profilepath STR] --runid STR --readtype
                  {2d,temp,comp,mixed,unknown} --readclass {all,pass,fail}
                  --infastqpath FILE [--targetrefpath FILE]
                  [--controlrefpath FILE] --outdir DIR --outprefix STR
                  [--useintdir BOOL] [--threads INT] --realmapprog
                  {bwa,graphmap} [--realmapparams STR] --randmapprog
                  {bwa,graphmap} [--randmapparams STR] [--overwrite BOOL]
                  [--savealignments BOOL] [--dryrun BOOL] [--touch BOOL]
                  [--ruffusgraph BOOL] [--verbose BOOL]
                  [--deleteintfiles BOOL] [--logindentwidth INT]
                  [--covplotwinsz INT] [--usecbcolours BOOL] [--prog_bwa STR]
                  [--prog_graphmap STR] [--prog_samtools STR]
                  [--prog_rscript STR]

Optional Arguments:

  -h, --help            show this help message and exit
  --bindir DIR          Absolute path of directory containing this program and
                        .ini and .profile config files (default: None)
  --profilepath STR     Path to BASH environment variables path (default:
                        None)
  --runid STR           RunID, which should not contain any spaces or unusual
                        punctuation (default: None)
  --readtype {2d,temp,comp,mixed,unknown}
                        Type of ONT reads (default: None)
  --readclass {all,pass,fail}
                        Class of ONT reads (default: None)
  --infastqpath FILE    Input reads path (FASTQ) (default: None)
  --targetrefpath FILE  Path to FASTA file for extracting reads of target
                        genome (default: None)
  --controlrefpath FILE
                        Reference control sample(s) (FASTA) (default: None)
  --outdir DIR          Output directory (default: None)
  --outprefix STR       Prefix for output files (default: None)
  --useintdir BOOL      Run pipeline to create output files in an intermediate
                        subdir first, then move them to the final directory if
                        the pipeline completes without error (default: False)
  --threads INT         Maximum number of threads to use (default: 1)
  --realmapprog {bwa,graphmap}
                        Mapping program to use for the real ONT target and/or
                        control reads: bwa, graphmap (default: bwa)
  --realmapparams STR   Additional parameters to pass to --realmapprog
                        surrounded by double-quotes if it contains spaces
                        (default: mem -x ont2d -M)
  --randmapprog {bwa,graphmap}
                        Mapping program to use for the randomised ONT control
                        reads: bwa, graphmap (default: graphmap)
  --randmapparams STR   Additional parameters to pass to --randmapprog
                        surrounded by double-quotes if it contains spaces
                        (default: -x nanopore -C)
  --overwrite BOOL      Print warning and overwrite any output files that may
                        already exist (default: False)
  --savealignments BOOL
                        Save each alignment to a separate FASTA file for
                        diagnostic purposes (default: False)
  --dryrun BOOL         Run without creating any output or changing any states
                        (default: False)
  --touch BOOL          Force Ruffus to touch all files only, making them seem
                        up to date (default: False)
  --ruffusgraph BOOL    Produce Ruffus dependency graph only, without running
                        pipeline (default: False)
  --verbose BOOL        Verbose diagnostic output (default: False)
  --deleteintfiles BOOL
                        Delete intermediate output files before exiting
                        (default: False)
  --logindentwidth INT  Print log messages from each new spawned process
                        indented by this many spaces (default: 4)
  --covplotwinsz INT    Print log messages from each new spawned process
                        indented by this many spaces (default: 100)
  --usecbcolours BOOL   Use colour-blind-friendly colours (default: False)
  --prog_bwa STR        BWA program (absolute pathname) (default: bwa)
  --prog_graphmap STR   GraphMap program (absolute pathname) (default:
                        graphmap)
  --prog_samtools STR   SAMtools program (absolute pathname) (default:
                        samtools)
  --prog_rscript STR    Rscript program (absolute pathname) (default: Rscript)

Examples:


Description:

A mapping-based pipeline to report:
  - genuine mappings of reads to a reference
  - quality statistics for reads
  - yield statistics for the run
  
Should work for data from:
  - single-chromosome, pure bacterial or viral samples
  - bacterial samples with multiple chromosome(s) and/or plasmid(s)
    where there may be regions of significant homology between parts
    of the chromosomes and plasmids
  - eukaryotic chromosomes
  - chimeric reads
  
Pipeline steps:
  1.  Create a single FASTA with one contig per chromosome/plasmid/spike-in
  2.  Align all reads to the set of target and control references using "bwa mem"
      and generate alignment statistics.
  3.  Let the subset of reads that map to the control reference(s) be
      part of the training set for class 'real' and generate alignment statistics.
  4.  Let the subset of reads that map to the control reference(s) be
      randomised in base order and quality to be the training set for class 'random'
      and generate alignment statistics.
  5.  Let the subset of reads that map to the target reference(s)
      be the set of alignments to be classified.
  6.  Classify the target alignments based on whether a linear regression
      of the principal components of some discriminating alignment statistics
      deem the alignment close to random alignments or non-random alignments.
  7.  Generate the final set of alignment, read and run statistics based on
      non-random alignments.
  8.  Save unique, chimeric and unmapped reads in BAM format (NOT IMPLEMENTED)
  9.  Save statistics report in HTML format (NOT IMPLEMENTED)
  
Input:
  - reads in FASTQ format
  - reference(s) for target sample(s) in a single FASTA file
  - reference(s) for control sample(s) in a single FASTA file
  
Output:
  - outdir/outprefix_alignstats.txt
  - outdir/outprefix_readstats.txt
  - outdir/outprefix_runstats.txt
  
Pre-conditions:
  - bindir must appear on the command-line to define where
    to find auxiliary programs and python packages.
  - profilepath is useful for storing BASH commands for changing the environment
    before running it SGE scheduling.
```
