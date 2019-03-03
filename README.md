
## Some scripts for PacBio sequencing data processing


### Publication
---

+ **Song WZ**, Thomas T*, Edwards R* (2018) Complete genome sequences of pooled genomic DNA from 10 marine bacteria using PacBio long-read sequencing, (unpublished)
+ Contact: Weizhi Song (songwz03@gmail.com)
+ Affiliation: The Centre for Marine Bio-Innovation (CMB), The University of New South Wales, Sydney, Australia


### Extract reads from SAM file
---

1. Help information

        python3 get_reads_from_sam.py -h

        arguments:
          -h, --help      show this help message and exit
          -sam            Input sam file
          -ctgs           Contig id list
          -option         Specify '1' to get reads mapped to provided contigs, or '0' to get reads not mapped to provided contigs
          -out            Output reads

1. Example commands below, an example of the ctg_ids.txt file was provided

        # get reads mapped to provided contigs
        $ python3 get_reads_from_sam.py -sam input.sam -ctg ctg_ids.txt -option 1 -out mapped_reads.fasta

        # get reads not mapped to provided contigs
        $ python3 get_reads_from_sam.py -sam input.sam -ctg ctg_ids.txt -option 0 -out unmapped_reads.fasta


### Purity assessment for diploid assemblies
---

1. Help information

        python3 get_purity.py -h

        arguments:
          -h, --help    show this help message and exit
          -r            folder holds reference sequences
          -x            extension of reference files
          -q            query sequences
          -n            number of reads to simulate
          -p            purity cut-off for reference assignment, default: 85
          -l            reads length, default: 250
          -i            insertion size, default: 500
          -m            the minimum number of mapped reads required for purity
                        calculation, NON-ZERO, default: 100
          -bbmap        path to BBMAP executable, default: bbmap.sh
          -circle       specify to simulate reads crossing the break point if the
                        references are circularized
          -split        export forward/reverse reads into separate files
          -keep_temp    keep temporary files
          -quiet        suppress reporting information


1. Example commands

        # assess the purity of query contigs with one million short reads simulated from provided reference genomes
        $ python3 get_purity.py -r ref_genome_folder -x fasta -q query_contigs.fasta -n 1000000

1. Output purity results

    |Sequence_ID|Assignment|Purity(%)|Length(Mbp)|
    |---|---|---|---|
    |Ctg_A2|Reference_1|100.0|0.06|
    |Ctg_A3|Reference_1|99.76|1.81|
    |Ctg_A4|Reference_2|100.0|0.05|
    |Ctg_A5|Reference_1|99.85|0.3|
    |Ctg_B2|Reference_2|100.0|0.06|
    |Ctg_B4|Reference_1|100.0|0.08|
    |Ctg_A1|Ambiguous|40.54/59.46|0.07|
    |Ctg_B1|Ambiguous|50.88/49.12|0.07|

1. Output plot

    ![purity_plot](images/DSM17395.haplotigs.purity.png)



