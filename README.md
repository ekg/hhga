# HHGA

![orb](https://raw.githubusercontent.com/ekg/hhga/master/images/orb.jpg)

## Have Haplotypes Genotypes Alleles

### A tool to prepare genome alignments and genotypes for consumption by the Vowpal Wabbit

[![Build Status](https://travis-ci.org/ekg/hhga.svg)](https://travis-ci.org/ekg/hhga)

## Overview

HHGA transforms genomic alignments and candidate haplotypes into a feature-space representation suitable for use by machine learning systems, in particular [Vowpal Wabbit](https://github.com/JohnLangford/vowpal_wabbit).
We take inspiration from genome visualizers that are often used to inspect candidate variant calls.
For example, `samtools tview` implements a text-based view of alignments in a particular region of the reference genome with encodings representing the read placement, direction, non-reference bases in the alignments, gaps (insertions and deletions are both shown faithfully using a gap character), and the qualities of mismatches. 


![tview](https://raw.githubusercontent.com/ekg/hhga/master/images/tview.png)

Additional tools provide mechanisms to convert these representations into human-readable formats (both text and image) for debugging.

HHGA can be thought of as an example decision synthesizer.
We use it in combination with genomic truth sets to generate examples from sequencing data that can be used in decision learning.
Our primary interest is to approximate P(genotype|alignments), where the genotype is composed of N haplotypes (typically 2), and the alignments represent the observational data we have at the genomic locus.

## Feature space representation

We represent the MSA of the reads and the reference, the estimated per-base errors, the read strands releative to the reference, the alignment mapping qualities, and the candidate haplotypes in the genotype of the sample at the locus using the namespaced libSVM format used by Vowpal Wabbit.

As an example, if this is a tview-like representation of a set of alignments at a putative SNP site, augmented with candidate haplotypes. Generated via `hhga -b minigiab/NA12878.chr22.tiny.bam -f minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -r q:10502-10562 -w 32 -t -x 5` from the `test/` directory.

```txt
q_10532_C_A
reference   AAAAACAAAAAAC-AACAAAAAAAAGGAAGGA
hap                         .               
hap                         A               
geno                        .               
geno                        A               
sOdqFxYPZI  .............-.................. 0.0	1.0 0.0 60 chr22.bin8.cram:166:8827
sOdqFxYPZI  .............-.................. 0.1	1.0 0.0 60 chr22.bin8.cram:166:8828
SodqFxYPZI  NNN..........-.................. 0.2	1.0 0.0 60 chr22.bin8.cram:166:8817
sOdqfXYPZI  .............-.................. 0.3	1.0 0.0 60 chr22.bin8.cram:166:8830
SodqFxYPZI  .............-.................. 0.4	1.0 0.0 60 chr22.bin8.cram:166:8777
SodqfXYPZI  .............-..A............... 1.0	0.0 1.0 60 chr22.bin8.cram:166:8789
SodqFxYPZI  .............-..A............... 1.1	0.0 1.0 60 chr22.bin8.cram:166:8794
sOdqfXYPZI  .............-..A........A...CA. 1.2	0.0 1.0 60 chr22.bin8.cram:166:8837
SodqFxYPZI  .............A..A............... 1.3	0.0 1.0 60 chr22.bin8.cram:166:8797
SodqfXYPZI  .............-..A............... 1.4	0.0 1.0 60 chr22.bin8.cram:166:8811
SodqFxYPZI  .............-..-............... 8.0	0.0 0.0 60 chr22.bin8.cram:166:8802
sOdqfXYPZI                   ............... 8.1	0.0 0.0 60 chr22.bin8.cram:166:8855
sOdqFxYPZI                         ......... 8.2	0.0 0.0 60 chr22.bin8.cram:166:8856
SodqfXYPZI  ........                         8.3	0.0 0.0 60 chr22.bin8.cram:166:8780
sOdqFxYPZI                            ...... 8.4	0.0 0.0 60 chr22.bin8.cram:166:8857
DPSum_1:705.0 HRun_1:2.0 HapNoVar_1:0.0 NoPLTot_1:0.0 PLminsumOverDP_1:12.6 PLminsum_1:8913.0 QUAL:8913.0 TrancheABQDmin2_1:0.0 TrancheAlignmin2_1:0.0 TrancheMapmin2_1:0.0 TrancheSSEmin2_1:0.0 YesPLtot_1:8.0 allalts_1:0.0 datasetcalls_1:8.0 genoMapGood_1:8.0 geno_1:2.0 platforms_1:2.0 PLCG_1:653 PLCG_2:0 PLCG_3:743 PLHSWG_1:1626 PLHSWG_2:0 PLHSWG_3:1610 PLILL250_1:501 PLILL250_2:0 PLILL250_3:624 PLILLCLIA_1:896 PLILLCLIA_2:0 PLILLCLIA_3:1250 PLILLWG_1:320 PLILLWG_2:0 PLILLWG_3:420 PLIllPCRFree_1:908 PLIllPCRFree_2:0 PLIllPCRFree_3:1146 PLPlatGen_1:3625 PLPlatGen_2:0 PLPlatGen_3:4125 PLXIll_1:384 PLXIll_2:0 PLXIll_3:581 platformbias_1:none platformnames_1:ill platformnames_2:cg varType_1:SNP 
```    

The feature space transformation of this site would be given by `hhga -b minigiab/NA12878.chr22.tiny.bam -f minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -r q:10502-10562 -w 32 -x 5``. The output will come on a single line.

For debugging we can pipe the hhga format output through `sed s/\|/\\n\|/g | column -t` to convert the spaces into newlines and columnarize the fields.

```txt
'q_10532_C_A
|ref            1A:1         2A:1      3A:1          4A:1         5A:1                    6C:1             7A:1       8A:1                 9A:1                  10A:1               11A:1               12A:1         13C:1        14U:1             15A:1            16A:1     17C:1          18A:1   19A:1   20A:1   21A:1   22A:1   23A:1   24A:1   25A:1   26G:1   27G:1   28A:1   29A:1   30G:1   31G:1   32A:1
|hap1           1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17R:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|hap2           1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17A:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|geno1          1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17R:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|geno2          1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17A:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln0.0         1R:40        2R:40     3R:40         4R:40        5R:40                   6R:40            7R:40      8R:40                9R:40                 10R:40              11R:40              12R:37        13R:27       14U:1             15R:37           16R:40    17R:40         18R:40  19R:37  20R:33  21R:37  22R:33  23R:40  24R:22  25R:40  26R:27  27R:37  28R:37  29R:40  30R:15  31R:37  32R:40
|aln0.1         1R:40        2R:40     3R:37         4R:40        5R:40                   6R:40            7R:40      8R:40                9R:40                 10R:40              11R:40              12R:40        13R:40       14U:1             15R:40           16R:40    17R:40         18R:40  19R:37  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:37  31R:40  32R:40
|aln0.2         1N:6         2N:6      3N:6          4R:40        5R:40                   6R:40            7R:40      8R:40                9R:40                 10R:40              11R:40              12R:40        13R:40       14U:1             15R:40           16R:40    17R:40         18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln0.3         1R:40        2R:40     3R:40         4R:33        5R:40                   6R:27            7R:40      8R:40                9R:37                 10R:22              11R:37              12R:40        13R:15       14U:1             15R:40           16R:40    17R:27         18R:37  19R:33  20R:33  21R:40  22R:40  23R:37  24R:37  25R:15  26R:15  27R:15  28R:37  29R:33  30R:27  31R:22  32R:37
|aln0.4         1R:40        2R:40     3R:40         4R:40        5R:40                   6R:27            7R:40      8R:40                9R:40                 10R:40              11R:37              12R:37        13R:15       14U:1             15R:40           16R:37    17R:22         18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:37  30R:40  31R:33  32R:27
|aln1.0         1R:40        2R:37     3R:37         4R:37        5R:15                   6R:22            7R:33      8R:40                9R:37                 10R:40              11R:40              12R:40        13R:15       14U:1             15R:40           16R:40    17A:40         18R:40  19R:40  20R:40  21R:37  22R:37  23R:40  24R:37  25R:40  26R:40  27R:33  28R:37  29R:40  30R:37  31R:37  32R:37
|aln1.1         1R:40        2R:40     3R:40         4R:40        5R:37                   6R:40            7R:40      8R:37                9R:40                 10R:40              11R:37              12R:40        13R:40       14U:1             15R:40           16R:33    17A:40         18R:40  19R:40  20R:40  21R:40  22R:37  23R:40  24R:40  25R:40  26R:40  27R:37  28R:15  29R:40  30R:40  31R:40  32R:37
|aln1.2         1R:33        2R:40     3R:40         4R:40        5R:40                   6R:33            7R:37      8R:27                9R:27                 10R:40              11R:40              12R:40        13R:15       14U:1             15R:40           16R:37    17A:40         18R:40  19R:40  20R:33  21R:33  22R:22  23R:40  24R:27  25R:15  26A:15  27R:15  28R:15  29R:37  30C:15  31A:15  32R:33
|aln1.3         1R:40        2R:40     3R:40         4R:40        5R:15                   6R:33            7R:40      8R:40                9R:37                 10R:40              11R:40              12R:33        13R:27       14A:27            15R:40           16R:40    17A:40         18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln1.4         1R:27        2R:27     3R:33         4R:37        5R:37                   6R:27            7R:40      8R:40                9R:40                 10R:40              11R:37              12R:37        13R:15       14U:1             15R:40           16R:40    17A:40         18R:40  19R:40  20R:37  21R:40  22R:40  23R:40  24R:37  25R:37  26R:40  27R:33  28R:37  29R:37  30R:40  31R:37  32R:15
|aln8.0         1R:40        2R:40     3R:37         4R:40        5R:37                   6R:15            7R:40      8R:40                9R:40                 10R:40              11R:40              12R:40        13R:33       14U:1             15R:40           16R:40    17U:40         18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:37  27R:40  28R:40  29R:37  30R:15  31R:40  32R:33
|aln8.1         1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17M:1          18R:33  19R:33  20R:37  21R:37  22R:37  23R:40  24R:40  25R:40  26R:33  27R:40  28R:37  29R:37  30R:37  31R:37  32R:40
|aln8.2         1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17M:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24R:33  25R:33  26R:27  27R:37  28R:37  29R:15  30R:33  31R:15  32R:40
|aln8.3         1R:40        2R:40     3R:37         4R:37        5R:15                   6R:37            7R:33      8R:33                9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17M:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln8.4         1M:1         2M:1      3M:1          4M:1         5M:1                    6M:1             7M:1       8M:1                 9M:1                  10M:1               11M:1               12M:1         13M:1        14M:1             15M:1            16M:1     17M:1          18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27R:37  28R:37  29R:37  30R:37  31R:37  32R:33
|col0           R:40         R:40      N:6           R:40         R:40                    R:40             R:40       R:33                 R:40                  R:27                R:40                M:1           M:1          R:40              M:1
|col1           R:40         R:40      N:6           R:40         R:40                    R:37             R:40       R:40                 R:40                  R:27                R:40                M:1           M:1          R:40              M:1
|col2           R:40         R:37      N:6           R:40         R:40                    R:37             R:40       R:40                 R:40                  R:33                R:37                M:1           M:1          R:37              M:1
|col3           R:40         R:40      R:40          R:33         R:40                    R:37             R:40       R:40                 R:40                  R:37                R:40                M:1           M:1          R:37              M:1
|col4           R:40         R:40      R:40          R:40         R:40                    R:15             R:37       R:40                 R:15                  R:37                R:37                M:1           M:1          R:15              M:1
|col5           R:40         R:40      R:40          R:27         R:27                    R:22             R:40       R:33                 R:33                  R:27                R:15                M:1           M:1          R:37              M:1
|col6           R:40         R:40      R:40          R:40         R:40                    R:33             R:40       R:37                 R:40                  R:40                R:40                M:1           M:1          R:33              M:1
|col7           R:40         R:40      R:40          R:40         R:40                    R:40             R:37       R:27                 R:40                  R:40                R:40                M:1           M:1          R:33              M:1
|col8           R:40         R:40      R:40          R:37         R:40                    R:37             R:40       R:27                 R:37                  R:40                R:40                M:1           M:1          M:1               M:1
|col9           R:40         R:40      R:40          R:22         R:40                    R:40             R:40       R:40                 R:40                  R:40                R:40                M:1           M:1          M:1               M:1
|col10          R:40         R:40      R:40          R:37         R:37                    R:40             R:37       R:40                 R:40                  R:37                R:40                M:1           M:1          M:1               M:1
|col11          R:37         R:40      R:40          R:40         R:37                    R:40             R:40       R:40                 R:33                  R:37                R:40                M:1           M:1          M:1               M:1
|col12          R:27         R:40      R:40          R:15         R:15                    R:15             R:40       R:15                 R:27                  R:15                R:33                M:1           M:1          M:1               M:1
|col13          U:1          U:1       U:1           U:1          U:1                     U:1              U:1        U:1                  A:27                  U:1                 U:1                 M:1           M:1          M:1               M:1
|col14          R:37         R:40      R:40          R:40         R:40                    R:40             R:40       R:40                 R:40                  R:40                R:40                M:1           M:1          M:1               M:1
|col15          R:40         R:40      R:40          R:40         R:37                    R:40             R:33       R:37                 R:40                  R:40                R:40                M:1           M:1          M:1               M:1
|col16          R:40         R:40      R:40          R:27         R:22                    A:40             A:40       A:40                 A:40                  A:40                U:40                M:1           M:1          M:1               M:1
|col17          R:40         R:40      R:40          R:37         R:40                    R:40             R:40       R:40                 R:40                  R:40                R:40                R:33          M:1          M:1               M:1
|col18          R:37         R:37      R:40          R:33         R:40                    R:40             R:40       R:40                 R:40                  R:40                R:40                R:33          M:1          M:1               M:1
|col19          R:33         R:40      R:40          R:33         R:40                    R:40             R:40       R:33                 R:40                  R:37                R:40                R:37          M:1          M:1               M:1
|col20          R:37         R:40      R:40          R:40         R:40                    R:37             R:40       R:33                 R:40                  R:40                R:40                R:37          M:1          M:1               M:1
|col21          R:33         R:40      R:40          R:40         R:40                    R:37             R:37       R:22                 R:40                  R:40                R:40                R:37          M:1          M:1               M:1
|col22          R:40         R:40      R:40          R:37         R:40                    R:40             R:40       R:40                 R:40                  R:40                R:40                R:40          M:1          M:1               M:1
|col23          R:22         R:40      R:40          R:37         R:40                    R:37             R:40       R:27                 R:40                  R:37                R:40                R:40          R:33         M:1               M:1
|col24          R:40         R:40      R:40          R:15         R:40                    R:40             R:40       R:15                 R:40                  R:37                R:40                R:40          R:33         M:1               M:1
|col25          R:27         R:40      R:40          R:15         R:40                    R:40             R:40       A:15                 R:40                  R:40                R:37                R:33          R:27         M:1               M:1
|col26          R:37         R:40      R:40          R:15         R:40                    R:33             R:37       R:15                 R:40                  R:33                R:40                R:40          R:37         M:1               R:37
|col27          R:37         R:40      R:40          R:37         R:40                    R:37             R:15       R:15                 R:40                  R:37                R:40                R:37          R:37         M:1               R:37
|col28          R:40         R:40      R:40          R:33         R:37                    R:40             R:40       R:37                 R:40                  R:37                R:37                R:37          R:15         M:1               R:37
|col29          R:15         R:37      R:40          R:27         R:40                    R:37             R:40       C:15                 R:40                  R:40                R:15                R:37          R:33         M:1               R:37
|col30          R:37         R:40      R:40          R:22         R:33                    R:37             R:40       A:15                 R:40                  R:37                R:40                R:37          R:15         M:1               R:37
|col31          R:40         R:40      R:40          R:37         R:27                    R:37             R:37       R:33                 R:40                  R:15                R:33                R:40          R:40         M:1               R:33
|match0.0       1H:1         2H:0
|match0.1       1H:1         2H:0
|match0.2       1H:1         2H:0
|match0.3       1H:1         2H:0
|match0.4       1H:1         2H:0
|match1.0       1H:0         2H:1
|match1.1       1H:0         2H:1
|match1.2       1H:0         2H:1
|match1.3       1H:0         2H:1
|match1.4       1H:0         2H:1
|match8.0       1H:0         2H:0
|match8.1       1H:0         2H:0
|match8.2       1H:0         2H:0
|match8.3       1H:0         2H:0
|match8.4       1H:0         2H:0
|properties0.0  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties0.1  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties0.2  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties0.3  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:0          xmate:1    ymap:1               paired:1              zprimary:1          iproper:1
|properties0.4  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties1.0  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:0          xmate:1    ymap:1               paired:1              zprimary:1          iproper:1
|properties1.1  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties1.2  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:0          xmate:1    ymap:1               paired:1              zprimary:1          iproper:1
|properties1.3  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties1.4  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:0          xmate:1    ymap:1               paired:1              zprimary:1          iproper:1
|properties8.0  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties8.1  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:0          xmate:1    ymap:1               paired:1              zprimary:1          iproper:1
|properties8.2  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|properties8.3  mapqual:60   strand:1  ostrand:0     dup:0        qcfail:0                fmate:0          xmate:1    ymap:1               paired:1              zprimary:1          iproper:1
|properties8.4  mapqual:60   strand:0  ostrand:1     dup:0        qcfail:0                fmate:1          xmate:0    ymap:1               paired:1              zprimary:1          iproper:1
|software       DPSum_1:705  HRun_1:2  HapNoVar_1:0  NoPLTot_1:0  PLminsumOverDP_1:12.64  PLminsum_1:8913  QUAL:8913  TrancheABQDmin2_1:0  TrancheAlignmin2_1:0  TrancheMapmin2_1:0  TrancheSSEmin2_1:0  YesPLtot_1:8  allalts_1:0  datasetcalls_1:8  genoMapGood_1:8  geno_1:2  platforms_1:2
```

We use six symbols to encode the MSA, `{ A, T, G, C, N, U, M, R, S }`, where:

* `A`, `T`, `G`, and `C` are DNA bases (by default only used in `|ref`)
* `N` is the degenerate symbol representing lack of information of actual base but knowledge of sequence length
* `U` is a gap symbol required to normalize the MSA into a matrix
* `M` is a symbol indicating if we don't have any information at the position in the MSA
* `R` to indicate when a base is the same as the reference base.
** This reduces the feature complexity of the model and helps performance. We further extend this to one specific symbol for each position in the matrix. So the reference-matching base at position 23 is always `23R`.
* `S` represents soft clips

`U`, and `M` always have weight 1. Reference bases always have weight 1. Bases in the alignments have weights based on their phred scores. The weight of `S` features is the length of the soft clip.

Various features of the reads are represented in other namespaces.
    
The first entry in the line defines the class of the example. The correspondence between the classes and genotypes is: 1: 0/0, 2: 0/1, 3: 0/2, 4:1/1, 5:1/2, 6:2/2, 7:unkown

* `ref` : the reference
* `hap*`: the haplotypes described in the VCF record (including reference allele)
* `geno*`: the genotypes described in the VCF record (e.g. we have 2N of these for dipliods)
* `0.0|, 0.1|, 1.0|, ... 9.0|` : the alignments overlapping the locus grouped by the haplotype(s) they match betst with with a namespace at the end for softclipped reads
* `match` : the rate of match between the alignment and haplotypes
* `properties` : things about the alignment taken from the input
* `software`: annotations from the VCF file, often software specific features

The reference sequence is given with the underlying bases in the `|ref` namespace.

Each alignment's namespaces in `|match` has these features:
    
* `*H` : the match between this alignment and the given haplotype from the VCF
    
Each alignment's namespaces in `|properties` has these features:

* `mapq` : p(mapping correct)
* `strand` : 1 if reversed, 0 otherwise (not written)
* `ostrand` : mate's strand, 1 if reversed
* `dup` : 1 if the read is marked as a duplicate
* `qcfail` : 1 if the read has failed some QC
* `fmate` : 1 if the read is the first mate
* `xmate` : 1 if the read is the second mate
* `ymap` : 1 if the mate is mapped (this mate must be mapped to be output)
* `paired` : 1 if the read is paired
* `zprimary` : 1 if the read is the "primary" alignment
* `iproper` : 1 if the read is in a "proper" pair, in the expected configuration and distance from its mate


## Usage

Sketch of usage. We extract the feature space representation using windows of size 100 at sites that are defined in the candidates in the variant input file. We use `truth.vcf.gz` to indicate which genotypes and alleles are true. All others in `vars.vcf.gz` are assumed false.

```bash
( hhga -v trues.vcf.gz  -w 100 -b aln.bam -f ref.fa -c 1 \
  hhga -v falses.vcf.gz -w 100 -b aln.bam -f ref.fa -c -1 ) \
    | shuf \
    | vw --save_resume -c --passes 7 --boosting 1 --ngram a3h3
```

## Building

```
git clone --recursive https://github.com/ekg/hhga.git
cd hhga
source ./source_me.sh
make
make test
```

The `hhga` executable is at `bin/hhga`.

## TODO

- add another symbol for the soft clipping
- add inputs beyond freebayes
- find how to express challenge objectives as vw parametrization
- 



## FDA Challenge objectives

"that all submitted accuracy comparisons reach a minimum threshold of 90% for each of the precision and the recall statistics."
Precision (PPV)	(true positives) / (true positives + false positives)
Recall (sensitivity)	(true positives) / (true positives + false negatives)
F-measure	2 * precision * recall / (precision + recall)

