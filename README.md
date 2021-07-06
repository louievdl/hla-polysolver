# louievdl/hla-polysolver

<a href="https://software.broadinstitute.org/cancer/cga/polysolver">Polysolver software</a> is developed by the Broad Institute. Availability and running instructions for the most current version are described <a href="https://software.broadinstitute.org/cancer/cga/polysolver_download">here</a>. This version is preferred.

This fork of Polysolver is a minor modification on <a href="https://github.com/jason-weirather/hla-polysolver">jason-weirather/hla-polysolver</a>, which 1) was based on Polysolver 1.0; 2) was modified under Polysolver's BSD-style License; and 3) does not produce exactly the same HLA haplotype calls as the example file included with Polysolver. Most differences were very minor, occurring in the 7th and 8th digits, with one difference in the 4th digit of the HLA-B allele. It should be noted that this work is not affiliated with Polysolver or the Broad Institute, and my work is free to use, reuse and modify under Apache License 2.0. See `LICENSE` file for licenses of software dependencies.

The purpose of this update of hla-polysolver is to: 1) handle any genome build (hg18/GRCh36, hg19/GRCh37, or hg38/GRCh38); 2) tolerate presence or absence of 'chr' in the chromosome specification of bam files; and 3) like <a href="https://github.com/jason-weirather/hla-polysolver">jason-weirather/hla-polysolver</a>, provide HLA haplotype "winner" files for loss of heterozygosity analysis by `LOHHLA` software (<a href="https://pubmed.ncbi.nlm.nih.gov/29107330">McGranahan <i>et al</i></a>). As such, it provides one part of a Polysolver-LOHHLA pipeline.

If you are using hla-polysolver in any form, you are requested to cite the original paper by <a href="https://www.ncbi.nlm.nih.gov/pubmed/26372948">Shukla <i>et al</i></a>.

## Software changes

#### polysolver v1.0 -> jason-weirather/hla-polysolver 1.0.0

* Added build recipe for conda
* Remove absolute path references that break the run
* Use the install of `build.sh` to make it so environment variables are set automatically **when run within Conda**
* Reduced use of environment variables
* Remove Novoalign index from the source code. Now is pre-built when building the conda environment, and is included along with other necessary data in the conda environment. If you are doing a local run and not using Conda see `build.sh` to see how to create this data file.
* Change from hardcoded perl and bash to whichever the user has installed under /usr/bin/env. This is better for non-conda runs but has the added benefit of using Conda's perl when running.
* Installed the shell scripts for hla typing, mutation calling, and annotation in the Conda enviornment so they are in the `PATH` of the Conda environment.
* Cleaning up the command calls and piping allows running the installed scripts from outside of the source directory.
* Removed hardcoded author paths
* Hardcoded temporary directory to /tmp. Not a great thing, but should work on most linux, and I plan to fix this soon.
* Added old picard tools dependency (likely what polysolver referred to as GATK)
* Updated data to include necessary fastas to complete mutation calling pipeline (part 2 of polysolver)

#### jason-weirather/hla-polysolver 1.0.0 -> louievdl/hla-polysolver 1.0.0

* allow user to specify any of the last three genome builds (hg18/GRCh36, hg19/GRCh37, or hg38/GRCh38)
* handle bam files with or without 'chr' in chromosome names
* removed hardwiring of HLA allele fasta file and index; user can provide an alternate. For consistency, the user can use the same allele file for other pipelines such as LOHHLA.
* HLA allele unique sequence file, allele IDs file, and SAM header file are all now generated on the fly from the allele fasta file

## A conda environment for use with hla-polysolver

Anaconda provides a package incorporating jason-weirather/hla-polsolver 1.0.0. Install and test it before overwriting with this fork. Disclaimer: as with any installation, success here may require supplying additional prerequisites specific to your system either at the level of the operating system or as conda packages. This is left as an exercise for the user.

Set up your conda channels in `~/.condarc`:
```
channels:
  - defaults
  - bioconda
  - conda-forge
  - vacation
```

Create a separate environment for polysolver to isolate its many old dependencies. Then activate and set the perl5 path.
```
$ conda create -n polysolver -c vacation hla-polysolver
$ conda activate polysolver
(polysolver)$ export PERL5LIB="$CONDA_PREFIX/lib/perl5/5.22.0/"
```

Get the github repository. You will access test bams and winner files now, and later use it to overwrite hla-polysolver in your conda installation
```
(polysolver)$ git clone https://github.com/louievdl/hla-polysolver.git
```

Run tests using anaconda's polysolver and bams from hla-polysolver repo. Increase java memory beforehand if necessary
The first script, `shell_call_hla_type` is the only strictly-necessary script for a LOHHLA pipeline, as it produces the winner files required by LOHHLA. It successfully processes the test files without modification.
The second script, `shell_call_hla_mutations_from_type`, to run successfully, requires minor modification, running muTect specifically under java version 7.
The third script, `shell_annotate_hla_mutations` runs on the test files without error.
```
(polysolver)$ export _JAVA_OPTIONS="-Xmx1g" && shell_call_hla_type hla-polysolver/test/test.bam Unknown 1 hg19 STDFQ 0 hla-polysolver/output
(polysolver)$ export _JAVA_OPTIONS="-Xmx1g" && shell_call_hla_mutations_from_type hla-polysolver/test/test.bam hla-polysolver/test/test.tumor.bam hla-polysolver/output/winners.hla.txt hg19 STDFQ hla-polysolver/output
(polysolver)$ shell_annotate_hla_mutations indiv hla-polysolver/output
```

Results will be in the hla-polysolver/output folder.

## HLA-POLYSOLVER MANUAL

#### TABLE OF CONTENTS ####
1. Description
    1.1 POLYSOLVER
    1.2 POLYSOLVER-based mutation detection
    1.3	Annotation of mutations
2. Installation
3. Testing
    3.1 POLYSOLVER
    3.2 POLYSOLVER-based mutation detection
    3.3	Annotation of mutations
4. Running
    4.1 POLYSOLVER
    4.2 POLYSOLVER-based mutation detection
    4.3	Annotation of mutations

#### 1. Description ####

This software package consists of 3 main tools:

1.1 POLYSOLVER (POLYmorphic loci reSOLVER)

This tool can be used for HLA typing based on an input exome BAM file and is currently infers infers alleles for the three major MHC class I  (HLA-A, -B, -C).

Script: shell_call_hla_type

Input parameters:
	
	-bam: path to the BAM file to be used for HLA typing
	-race: ethnicity of the individual (Caucasian, Black, Asian or Unknown)
	-includeFreq: flag indicating whether population-level allele frequencies should be used as priors (0 or 1)
	-build: reference genome used in the BAM file (hg18 or hg19)
	-format: fastq format (STDFQ, ILMFQ, ILM1.8 or SLXFQ; see Novoalign documentation)
	-insertCalc: flag indicating whether empirical insert size distribution should be used in the model (0 or 1)
	-outDir: output directory

Output:

	winners.hla.txt: file containing the two inferred alleles for each of HLA-A, HLA-B and HLA-C
	  
	  
1.2 POLYSOLVER-based mutation detection

This tool works on a tumor/normal pair of exome BAM files and inferred mutations in the tumor file. It assumes that POLYSOLVER has already been run on the normal BAM.

Script: shell_call_hla_mutations_from_type

Input parameters:

	-normal_bam_hla: path to the normal BAM file
	-tumor_bam_hla: path to the tumor BAM file
	-hla: inferred HLA allele file from POLYSOLVER (winners.hla.txt or winners.hla.nofreq.txt)
	-build: reference genome used in the BAM file (hg18 or hg19)
	-format: fastq format (STDFQ, ILMFQ, ILM1.8 or SLXFQ; see Novoalign documentation)
	-outDir: output directory	  
	
Output:

	call_stats.$allele.out: Mutect output for each inferred allele in winners.hla.txt
	$allele.all.somatic.indels.vcf: Strelka output for each inferred allele in winners.hla.txt
	
1.3 Annotation of mutations

This tool annotates the predicted mutations from (ii) with gene compartment and amino acid change information

Script: shell_annotate_hla_mutations

Input parameters:

	-indiv: individual ID, used as prefix for output files
	-dir: directory containing the raw call files (Mutect: call_stats*, Strelka: *all.somatic.indels.vcf). Also the output directory	
	
Output:

(a). Mutect
	$indiv.mutect.unfiltered.nonsyn.annotated -  list of all unfiltered mutations
	$indiv.mutect.filtered.nonsyn.annotated -  list of cleaned non-synonymous mutations
	$indiv.mutect.filtered.syn.annotated - list of cleaned synonymous changes
	$indiv.mutect.ambiguous.annotated - list of ambiguous calls. This will generally be empty (save for the header). It will be populated if the same mutation (ex. p.A319E) is found in two or more alleles in the individual, with the same allele fractions. In such cases one allele is randomly chosen and included in the .nonysn.annotated file while the complete list of alleles is listed in the .ambiguous.annotated file. If the ethnicity of the individual is known, an alternate method would be to pick the allele with the highest frequency.

(b). Strelka
	$indiv.mutect.unfiltered.nonsyn.annotated -  list of all unfiltered indels (as detected by Strelka)
	$indiv.strelka_indels.filtered.annotated - list of cleaned indels (as detected by Strelka)
	$indiv.strelka_indels.ambiguous.annotated - see description of $indiv.mutect.ambiguous.annotated in (a). above
	
	
#### 2. Installation ####

The POLYSOLVER suite of tools depends upon the following packages and utilities:


Samtools (http://samtools.sourceforge.net/)
GATK (https://www.broadinstitute.org/gatk/download)
Novoalign (http://www.novocraft.com/main/downloadpage.php)
Perl modules ((http://www.cpan.org/modules/INSTALL.html)
 - Math::BaseCalc
 - List::MoreUtils
 - List::Util
 - Parallel::ForkManager
 - POSIX
 - Dumpvalue
 - Data::Dumper
Bioperl (http://www.bioperl.org/wiki/Installing_BioPerl)
Mutect (http://www.broadinstitute.org/cancer/cga/mutect_download)
Strelka (https://sites.google.com/site/strelkasomaticvariantcaller/home/download)
 
Also make changes to the config.sh file to set up the following environmental variables

 -PSHOME: POLYSOLVER home directory
 -SAMTOOLS_DIR: directory containing the samtools executable
 -JAVA_DIR: directory containing the JAVA executable
 -NOVOALIGN_DIR: directory containing the Novoalign executables
 -GATK_DIR: directory containing the GATK jar files
 -MUTECT_DIR: directory containing the Mutect executable (for POLYSOLVER-based mutation detection only)
 -STRELKA_DIR: directory containing the Strelka  (for POLYSOLVER-based mutation detection only)

The following command should make the necessary changes prior to running the tools (assuming the tcsh shell):

source scripts/config.sh
  
 
#### 3. Testing ####

Your installation can be tested by running the following command from $PSHOME:

3.1 POLYSOLVER

scripts/shell_call_hla_type test/test.bam Unknown 1 hg19 STDFQ 0 test
 
If successful, the following command should not yield any differences:
 
diff test/winners.hla.txt test/orig.winners.hla.txt

3.2 POLYSOLVER-based mutation detection

scripts/shell_call_hla_mutations_from_type test/test.bam test/test.tumor.bam test/winners.hla.txt hg19 STDFQ test

If successful, the following command should not yield any differences:
 
diff test/call_stats.hla_b_39_01_01_02l.out test/orig.call_stats.hla_b_39_01_01_02l.out 

3.3 Annotation of mutations

scripts/shell_annotate_hla_mutations indiv test

If successful, the following command should not yield any differences:

diff test/indiv.mutect.filtered.nonsyn.annotated test/orig.indiv.mutect.filtered.nonsyn.annotated

#### 4. Running #####

The tools can be run using the following commands:

4.1 POLYSOLVER

$PSHOME/scripts/shell_call_hla_type </path/to/bam> <race> <includeFreq> <build> <format> <insertCalc> </path/to/output_directory>

example:

$PSHOME/scripts/shell_call_hla_type test/test.bam Unknown 1 hg19 STDFQ 0 test

4.2 POLYSOLVER-based mutation detection

$PSHOME/scripts/shell_call_hla_mutations_from_type </path/to/normal_bam> </path/to/tumor_bam> </path/to/winners.hla.txt> <build> <format> </path/to/output_directory>

example:

$PSHOME/scripts/shell_call_hla_mutations_from_type test/test.bam test/test.tumor.bam test/winners.hla.txt hg19 STDFQ test
 
4.3 Annotation of mutations

$PSHOME/scripts/shell_annotate_hla_mutations <prefix_to_use> </path/to/directory_with_mutation_detection_output>

example:

$PSHOME/scripts/shell_annotate_hla_mutations indiv test


