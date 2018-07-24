calcRoiCovg
===========
Given two BAM files, a reference sequence, and regions of interest, count the user customized
bp class types, for example, "AT,CpG,CG" sites with sufficient read-depth in both BAMs.

Usage
-----

Version 0.1
Usage: calcRoiCovg <bam1> <bam2> <roi_file> <ref_seq_fasta> <output_file>

        -q INT    filtering reads with mapping quality less than INT [20]
        -n INT    minimum reads depth for bam1 [6]
        -t INT    minimum reads depth for bam2 [8]
        -c STRING bp class types, delimited by comma, default: "AT,CpG,GC"


ROI file should be a tab-delimited list of [chrom, start, stop, annotation]
where start and stop are both 1-based chromosomal loci. For example:

20      44429404        44429608        ELMO2
MT      5903    7445    MT-CO1

ROI file *must* be sorted by chromosome/contig names


This tool was originally designed to count base-pairs that have sufficient read-depth for variant
calling across two BAM files (case vs control). The base-pairs are further classified into AT, CG
(non-CpG), and CpG sites, or whatever user custermized bp types with respect to the provided reference sequence. 
The resulting coverage stats are reported for each region of interest (ROI).

Sample ROI files can be found under the 'data' subdirectory. Note that they use 1-based loci, and
*must* be sorted by chromosome or contig names.

Install
-------
The Makefile assumes that you have the samtools source code in an environment variable `$SAMDIR`. If
you don't know what that means, then simply follow these steps from any directory that you have
permissions to write into:

Install some prerequisite packages if you are using Debian or Ubuntu:

    sudo apt-get install git libbam-dev zlib1g-dev

If you are using Fedora, CentOS or RHEL, you'll need these packages instead:

    sudo yum install git samtools-devel zlib-devel

Download the samtools-0.1.19 from SOURCEFORGE (http://sourceforge.net/projects/samtools/files/samtools/0.1.19):

    tar jxf samtools-0.1.19.tar.bz2
    cd samtools-0.1.19
    make
    export SAMTOOLS_ROOT=$PWD

Clone the calc-roi-covg repo, and build the `calcRoiCovg` binary:

    git clone https://github.com/Beifang/calcRoiCovg.git
    cd calc-roi-covg
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv calcRoiCovg /usr/local/bin/

xxx

