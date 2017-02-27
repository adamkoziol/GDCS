# geneSipprV2
==============
#Introduction

This pipeline searches for gene targets in FASTQ files. These files may be previously generated, or created from BCL 
files from a run of an Illumina MiSeq as part of the pipeline. The latter functionality is to allow for the creation
of FASTQ files from an in progress MiSeq run.


#Contents
This pipeline includes a main script (geneSipprV2, and several helper modules located in helper scripts, including

* fastqCreator
    * Creates FASTQ files from a provided Illumina MiSeq folder

* 16S
    * Performs genus level typing of prokaryotic samples

* MLST
    * Using the genus determined by 16S typing, the appropriate MLST scheme is used on the strain of interest. Currently
    only a select genera are represented, but it is possible to add any scheme desired

* pathotyping
    * CGE pathotype database + CHAS probes

* serotyping
    * CGE serotype database (currently only for *Escherichia* and *Salmonella*)

* virulence typing
    * CGE virulence gene database

* antimicrobial resistance marker finding
    * Uses ARMI

* rMLST
    * ribosomal multilocus sequence typing

* custom target analysis
    * Analysis of custom targets

#Arguments
-p, --path, Required. Path in which to place the reports folder. It is also used to find sequences and targets folders if arguments for these folders are not provided

-s, --sequencepath, Path of .fastq(.gz) files to process. If not provided, the default path of "path/sequences" will be used

-t, --targetpath, Path of target files to process. If not provided, the default path of "path/targets" will be used

-m, --miSeqPath, Path of the folder containing MiSeq run data folder (e.g. /media/miseq/MiSeqOutput)

-f, --miseqfolder, Name of the folder containing MiSeq run data (e.g. 151218_M02466_0126_000000000-AKF4P)

-r1, --readLengthForward, Length of forward reads to use. Can specify "full" to take the full length of forward reads specified on the SampleSheet

-r2, --readLengthReverse, Length of reverse reads to use. Can specify "full" to take the full length of reverse reads specified on the SampleSheet. Default is 0

-c, --customSampleSheet, Path of folder containing a custom sample sheet (still must be named "SampleSheet.csv)"

-P, --projectName, A name for the analyses. If nothing is provided, then the "Sample_Project" field in the provided sample sheet will be used. Please note that bcl2fastq 
     creates subfolders using the project name, so if multiple names are provided, the results will be split as into multiple projects
     
-16S, --16Styping, Perform 16S typing. Note that for analyses such as MLST, pathotyping, serotyping, and virulence typing that require the genus of a strain to proceed, 
     16S typing will still be performed
     
-M, --Mlst, Perform MLST analyses

-Y, --pathotYping, Perform pathotyping analyses

-S, --Serotyping, Perform serotyping analyses

-V, --Virulencetyping, Perform virulence typing analyses

-a, --armi, Perform ARMI antimicrobial typing analyses

-r, --rmlst, Perform rMLST analyses

-d, --detailedReports, Provide detailed reports with percent identity and depth of coverage values rather than just "+" for positive results

-C, --customTargetPath, Provide the path for a folder of custom targets .fasta format. Does not support multifasta files

#Use
Run geneSipprV2.py from the console with the desired arguments.

#Requirements
* Linux
* Python
* BioPython
* pysamstats and dependencies
* SMALT v?
* samtools
* mirabait
* more coming soon!

# Outputs
This pipeline generates multiple reports. They should be in the path/reports folder. A new folder is created for each analysis with the format: YYYY.MM.DD.HH.MM.SS