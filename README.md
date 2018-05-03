
Some scripts used for PacBio sequencing data processing
---

Publication
---
+ In preparation
+ Contact: Weizhi Song (songwz03@gmail.com)
+ Affiliation: The Centre for Marine Bio-Innovation (CMB), The University of New South Wales, Sydney, Australia


Extract needed reads from SAM file
---
        python3 get_reads_from_sam.py -h

        arguments:
          -h, --help      show this help message and exit
          -sam            Input sam file
          -ctgs           Contig id list
          -option         Specify '1' to get reads mapped to provided contigs, or '0' to get unmapped reads
          -out            Output reads


Purity assessment for diploid assemblies
---

        python3 get_purity.py -h

        arguments:
          -h, --help    show this help message and exit
          -r R          folder holds reference sequences
          -x X          the extension of reference files
          -q Q          query sequences
          -n N          the total number of reads need to be simulated
          -p P          purity cut-off for reference assignment, default: 85
          -l L          the length of simulated reads, default: 250
          -i I          the insert size of simulated reads, default: 500
          -m M          the minimum number of mapped reads required for purity
                        calculation, NON-ZERO, default: 500
          -bbmap BBMAP  path to BBMAP executable file, default: bbmap.sh
          -circle       specify to simulate reads crossing the break point if the
                        reference sequences are circularized
          -split        specify to export simulated forward and reverse reads into
                        separate files
          -keep_temp    specify to keep temporary files
          -quiet        specify to suppress reporting information

Output files
---

1. Purity plot

![purity_plot](images/DSM17395.haplotigs.purity.png)

1. Purity and reference assignment for each query sequences

        Sequence_ID	Assignment	Purity(%)	Length(Mbp)
        Ctg_A2	Reference_1	100.0	0.06
        Ctg_A3	Reference_1	99.76	1.81
        Ctg_A4	Reference_2	100.0	0.05
        Ctg_A5	Reference_1	99.85	0.3
        Ctg_B2	Reference_2	100.0	0.06
        Ctg_B4	Reference_1	100.0	0.08
        Ctg_A1	Ambiguous	40.54/59.46	0.07
        Ctg_B1	Ambiguous	50.88/49.12	0.07

