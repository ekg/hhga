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

As an example, if this is a tview-like representation of a set of alignments at a putative SNP site, augmented with candidate haplotypes. Generated via `hhga -b minigiab/9251-9252.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -w 32 -t` from the `test/` directory.

```txt
q_9251_GTTCT_G
reference   AGAAAGATTGTGCCAAGTTCTTTCTTTTTCAG
hap                         .....           
hap                         .----           
geno                        .....           
geno                        .----           
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8398
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8321
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8350
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8354
sOdqfXYPZI  .................----...G....... 0.2 1.0 60 chr22.bin8.cram:166:8443
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8392
sOdqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8446
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8447
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8448
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8449
sOdqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8450
sOdqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8451
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8393
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8448
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8394
sOdqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8455
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8390
sOdqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8457
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8436
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8400
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8460
sOdqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8461
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8462
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8428
sOdqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8464
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8420
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8426
SodqfXYPZI  ........................A..C.... 1.0 0.2 60 chr22.bin8.cram:166:8364
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8379
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8461
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8415
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8372
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8407
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8421
SodqFxYPZI  ...............................  1.0 0.2 60 chr22.bin8.cram:166:8377
sOdqfXYPZI     ............................. 1.0 0.2 60 chr22.bin8.cram:166:8474
SodqfXYPZI  ..........................       1.0 0.2 60 chr22.bin8.cram:166:8408
sOdqfXYPZI  ........................         1.0 0.2 60 chr22.bin8.cram:166:8436
SodqfXYPZI  .....................            1.0 0.2 60 chr22.bin8.cram:166:8361
sOdqfXYPZI             ..................... 1.0 0.2 60 chr22.bin8.cram:166:8475
SodqFxYPZI  ....................             0.8 0.2 60 chr22.bin8.cram:166:8341
sOdqFxYPZI  ...................              0.6 0.2 60 chr22.bin8.cram:166:8434
SodqFxYPZI                      G........... 0.0 0.0 60 chr22.bin8.cram:166:8432
SodqfXYPZI                      ............ 0.2 0.0 60 chr22.bin8.cram:166:8427
sOdqfXYPZI                      S........... 0.0 0.0 60 chr22.bin8.cram:166:8478
sOdqFxYPZI                      S........... 0.0 0.0 60 chr22.bin8.cram:166:8479
SodqFxYPZI                      S........... 0.0 0.0 60 chr22.bin8.cram:166:8443
SodqfXYPZI  ...........                      0.0 0.0 60 chr22.bin8.cram:166:8410
SodqfXYPZI  ...........                      0.0 0.0 60 chr22.bin8.cram:166:8409
SodqFxYPZI                         ......... 0.0 0.0 60 chr22.bin8.cram:166:8395
SodqfXYPZI  ........                         0.0 0.0 60 chr22.bin8.cram:166:8383
sOdqFxYPZI                          ........ 0.0 0.0 60 chr22.bin8.cram:166:8482
sOdqFxYPZI  .......                          0.0 0.0 60 chr22.bin8.cram:166:8428
sOdqFxYPZI  .                                0.0 0.0 60 chr22.bin8.cram:166:8427
SodqFxYPZI                                 . 0.0 0.0 60 chr22.bin8.cram:166:8424
AC_1:1.0 QUAL:359.9 
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

The feature space transformation of this site would be given by `hhga -b minigiab/9251-9252.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -w 20 -c 1`. The output will come on a single line.

```txt
q_9251_GTTCT_G
reference   AGAAAGATTGTGCCAAGTTCTTTCTTTTTCAG
hap                         .....           
hap                         .----           
geno                        .....           
geno                        .----           
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8398
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8321
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8350
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8354
sOdqfXYPZI  .................----...G....... 0.2 1.0 60 chr22.bin8.cram:166:8443
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8392
sOdqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8446
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8447
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8448
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8449
sOdqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8450
sOdqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8451
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8393
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8448
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8394
sOdqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8455
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8390
sOdqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8457
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8436
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8400
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8460
sOdqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8461
sOdqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8462
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8428
sOdqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8464
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8420
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8426
SodqfXYPZI  ........................A..C.... 1.0 0.2 60 chr22.bin8.cram:166:8364
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8379
SodqFxYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8461
SodqFxYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8415
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8372
SodqfXYPZI  ................................ 1.0 0.2 60 chr22.bin8.cram:166:8407
SodqfXYPZI  .................----........... 0.2 1.0 60 chr22.bin8.cram:166:8421
SodqFxYPZI  ...............................  1.0 0.2 60 chr22.bin8.cram:166:8377
sOdqfXYPZI     ............................. 1.0 0.2 60 chr22.bin8.cram:166:8474
SodqfXYPZI  ..........................       1.0 0.2 60 chr22.bin8.cram:166:8408
sOdqfXYPZI  ........................         1.0 0.2 60 chr22.bin8.cram:166:8436
SodqfXYPZI  .....................            1.0 0.2 60 chr22.bin8.cram:166:8361
sOdqfXYPZI             ..................... 1.0 0.2 60 chr22.bin8.cram:166:8475
SodqFxYPZI  ....................             0.8 0.2 60 chr22.bin8.cram:166:8341
sOdqFxYPZI  ...................              0.6 0.2 60 chr22.bin8.cram:166:8434
SodqFxYPZI                      G........... 0.0 0.0 60 chr22.bin8.cram:166:8432
SodqfXYPZI                      ............ 0.2 0.0 60 chr22.bin8.cram:166:8427
sOdqfXYPZI                      S........... 0.0 0.0 60 chr22.bin8.cram:166:8478
sOdqFxYPZI                      S........... 0.0 0.0 60 chr22.bin8.cram:166:8479
SodqFxYPZI                      S........... 0.0 0.0 60 chr22.bin8.cram:166:8443
SodqfXYPZI  ...........                      0.0 0.0 60 chr22.bin8.cram:166:8410
SodqfXYPZI  ...........                      0.0 0.0 60 chr22.bin8.cram:166:8409
SodqFxYPZI                         ......... 0.0 0.0 60 chr22.bin8.cram:166:8395
SodqfXYPZI  ........                         0.0 0.0 60 chr22.bin8.cram:166:8383
sOdqFxYPZI                          ........ 0.0 0.0 60 chr22.bin8.cram:166:8482
sOdqFxYPZI  .......                          0.0 0.0 60 chr22.bin8.cram:166:8428
sOdqFxYPZI  .                                0.0 0.0 60 chr22.bin8.cram:166:8427
SodqFxYPZI                                 . 0.0 0.0 60 chr22.bin8.cram:166:8424
AC_1:1.0 QUAL:359.9 

```

For debugging we can pipe the hhga format output through `sed s/\|/\\n\|/g | column -t` to convert the spaces into newlines and columnarize the fields.

```txt
1              'q_9251_GTTCT_G
|ref           1A:1             2G:1          3A:1       4A:1   5A:1      6G:1     7A:1     8T:1    9T:1      10G:1       11T:1      12G:1   13C:1   14C:1   15A:1   16A:1   17G:1   18T:1   19T:1   20C:1   21T:1   22T:1   23T:1   24C:1   25T:1   26T:1   27T:1   28T:1   29T:1   30C:1   31A:1   32G:1
|hap1          1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17R:1   18R:1   19R:1   20R:1   21R:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|hap2          1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17R:1   18U:1   19U:1   20U:1   21U:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|geno1         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17R:1   18R:1   19R:1   20R:1   21R:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|geno2         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17R:1   18U:1   19U:1   20U:1   21U:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln0          1R:37            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:37     10R:40      11R:37     12R:22  13R:37  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:40  21R:40  22R:37  23R:37  24R:37  25R:40  26R:40  27R:37  28R:37  29R:37  30R:37  31R:37  32R:37
|aln1          1R:40            2R:37         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:15  14R:37  15R:40  16R:40  17R:40  18R:40  19R:37  20R:40  21R:37  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:37  30R:37  31R:37  32R:37
|aln2          1R:40            2R:33         3R:37      4R:40  5R:40     6R:40    7R:37    8R:40   9R:40     10R:37      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18R:37  19R:40  20R:37  21R:40  22R:40  23R:40  24R:40  25R:37  26R:27  27R:27  28R:40  29R:40  30R:40  31R:22  32R:40
|aln3          1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:33  25R:40  26R:40  27R:40  28R:40  29R:37  30R:40  31R:40  32R:33
|aln4          1R:22            2R:15         3R:37      4R:33  5R:22     6R:40    7R:33    8R:27   9R:33     10R:22      11R:33     12R:33  13R:37  14R:37  15R:37  16R:27  17R:40  18U:37  19U:27  20U:40  21U:33  22R:33  23R:33  24R:37  25G:15  26R:37  27R:40  28R:22  29R:33  30R:27  31R:15  32R:22
|aln5          1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:37      11R:40     12R:40  13R:40  14R:40  15R:37  16R:33  17R:33  18U:37  19U:33  20U:33  21U:40  22R:40  23R:37  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:33  31R:40  32R:40
|aln6          1R:40            2R:40         3R:37      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18U:40  19U:40  20U:40  21U:40  22R:40  23R:40  24R:37  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln7          1R:37            2R:40         3R:37      4R:27  5R:40     6R:37    7R:40    8R:40   9R:40     10R:37      11R:37     12R:37  13R:40  14R:40  15R:37  16R:37  17R:40  18U:37  19U:37  20U:40  21U:40  22R:40  23R:37  24R:40  25R:40  26R:37  27R:40  28R:40  29R:40  30R:33  31R:40  32R:40
|aln8          1R:37            2R:40         3R:40      4R:40  5R:33     6R:37    7R:40    8R:33   9R:40     10R:40      11R:37     12R:40  13R:40  14R:33  15R:37  16R:40  17R:27  18U:37  19U:40  20U:27  21U:37  22R:37  23R:37  24R:33  25R:40  26R:40  27R:40  28R:40  29R:27  30R:40  31R:40  32R:37
|aln9          1R:40            2R:40         3R:37      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18U:40  19U:40  20U:40  21U:40  22R:40  23R:37  24R:40  25R:37  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln10         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:37    8R:37   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18U:40  19U:40  20U:40  21U:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln11         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:37  13R:40  14R:33  15R:40  16R:40  17R:40  18R:40  19R:37  20R:40  21R:37  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln12         1R:40            2R:37         3R:40      4R:40  5R:22     6R:40    7R:37    8R:40   9R:37     10R:40      11R:27     12R:37  13R:40  14R:40  15R:40  16R:40  17R:40  18U:40  19U:40  20U:40  21U:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:37  30R:40  31R:40  32R:40
|aln13         1R:40            2R:37         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:33  17R:40  18U:40  19U:33  20U:40  21U:37  22R:37  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln14         1R:33            2R:40         3R:37      4R:22  5R:33     6R:37    7R:37    8R:40   9R:33     10R:40      11R:40     12R:37  13R:27  14R:40  15R:37  16R:40  17R:37  18U:37  19U:40  20U:37  21U:37  22R:37  23R:33  24R:40  25R:27  26R:40  27R:37  28R:33  29R:37  30R:33  31R:40  32R:15
|aln15         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:33  17R:40  18U:40  19U:33  20U:40  21U:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln16         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:27   9R:37     10R:40      11R:40     12R:40  13R:40  14R:40  15R:37  16R:40  17R:40  18U:37  19U:40  20U:40  21U:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:37  30R:40  31R:40  32R:37
|aln17         1R:40            2R:15         3R:15      4R:40  5R:40     6R:37    7R:40    8R:40   9R:33     10R:33      11R:40     12R:40  13R:40  14R:40  15R:15  16R:33  17R:27  18R:33  19R:15  20R:37  21R:40  22R:40  23R:33  24R:40  25R:22  26R:40  27R:27  28R:40  29R:40  30R:40  31R:40  32R:40
|aln18         1R:40            2R:37         3R:40      4R:40  5R:40     6R:40    7R:37    8R:33   9R:40     10R:37      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:37  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln19         1R:40            2R:37         3R:40      4R:33  5R:37     6R:40    7R:40    8R:40   9R:37     10R:37      11R:40     12R:37  13R:40  14R:37  15R:37  16R:40  17R:37  18U:37  19U:40  20U:37  21U:40  22R:40  23R:37  24R:40  25R:37  26R:40  27R:40  28R:40  29R:40  30R:40  31R:37  32R:40
|aln20         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18U:40  19U:40  20U:40  21U:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln21         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln22         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18U:40  19U:40  20U:40  21U:37  22R:37  23R:40  24R:40  25R:37  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:37
|aln23         1R:37            2R:33         3R:37      4R:37  5R:33     6R:40    7R:33    8R:37   9R:40     10R:40      11R:37     12R:40  13R:40  14R:33  15R:40  16R:33  17R:40  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:37  26R:40  27R:37  28R:37  29R:40  30R:37  31R:37  32R:37
|aln24         1R:40            2R:37         3R:37      4R:40  5R:40     6R:37    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:37  21R:37  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:37  32R:40
|aln25         1R:22            2R:37         3R:40      4R:40  5R:40     6R:33    7R:40    8R:37   9R:40     10R:40      11R:40     12R:33  13R:40  14R:40  15R:40  16R:37  17R:33  18U:40  19U:37  20U:33  21U:40  22R:40  23R:40  24R:22  25R:27  26R:40  27R:40  28R:40  29R:40  30R:40  31R:37  32R:37
|aln26         1R:37            2R:37         3R:37      4R:37  5R:40     6R:27    7R:33    8R:15   9R:40     10R:37      11R:37     12R:37  13R:40  14R:40  15R:27  16R:40  17R:33  18R:15  19R:40  20R:27  21R:40  22R:27  23R:40  24R:40  25R:15  26R:33  27R:40  28R:40  29R:37  30R:40  31R:33  32R:40
|aln27         1R:15            2R:15         3R:40      4R:15  5R:37     6R:37    7R:27    8R:15   9R:15     10R:22      11R:15     12R:33  13R:22  14R:15  15R:40  16R:27  17R:15  18R:15  19R:27  20R:27  21R:15  22R:15  23R:15  24R:27  25A:15  26R:15  27R:22  28C:15  29R:15  30R:15  31R:27  32R:15
|aln28         1R:37            2R:15         3R:37      4R:37  5R:37     6R:37    7R:33    8R:37   9R:40     10R:33      11R:37     12R:22  13R:33  14R:37  15R:37  16R:27  17R:27  18R:22  19R:40  20R:27  21R:40  22R:33  23R:37  24R:40  25R:40  26R:27  27R:37  28R:33  29R:40  30R:37  31R:15  32R:33
|aln29         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln30         1R:37            2R:37         3R:37      4R:37  5R:40     6R:37    7R:37    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:37  16R:37  17R:33  18U:37  19U:37  20U:33  21U:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:33  31R:40  32R:27
|aln31         1R:37            2R:40         3R:27      4R:37  5R:33     6R:22    7R:37    8R:40   9R:27     10R:37      11R:37     12R:33  13R:40  14R:33  15R:37  16R:37  17R:40  18U:37  19U:37  20U:40  21U:37  22R:37  23R:33  24R:15  25R:37  26R:40  27R:27  28R:40  29R:37  30R:33  31R:40  32R:37
|aln32         1R:40            2R:37         3R:40      4R:37  5R:27     6R:40    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:40  17R:33  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:37  32R:40
|aln33         1R:40            2R:33         3R:37      4R:37  5R:27     6R:15    7R:27    8R:33   9R:33     10R:22      11R:33     12R:27  13R:37  14R:37  15R:15  16R:37  17R:37  18U:15  19U:37  20U:37  21U:33  22R:33  23R:37  24R:37  25R:33  26R:40  27R:40  28R:40  29R:33  30R:37  31R:37  32R:37
|aln34         1R:33            2R:40         3R:40      4R:37  5R:37     6R:40    7R:40    8R:37   9R:40     10R:40      11R:40     12R:40  13R:33  14R:40  15R:40  16R:37  17R:40  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:37  26R:40  27R:37  28R:37  29R:37  30R:33  31R:33  32M:1
|aln35         1M:1             2M:1          3M:1       4R:33  5R:37     6R:37    7R:37    8R:37   9R:40     10R:40      11R:40     12R:40  13R:33  14R:40  15R:40  16R:40  17R:40  18R:40  19R:40  20R:37  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln36         1R:33            2R:40         3R:33      4R:40  5R:40     6R:40    7R:37    8R:40   9R:40     10R:37      11R:37     12R:37  13R:37  14R:37  15R:40  16R:27  17R:40  18R:37  19R:40  20R:40  21R:22  22R:37  23R:37  24R:37  25R:33  26R:27  27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln37         1R:27            2R:22         3R:40      4R:33  5R:33     6R:33    7R:40    8R:40   9R:22     10R:37      11R:37     12R:37  13R:37  14R:37  15R:40  16R:37  17R:37  18R:27  19R:33  20R:40  21R:40  22R:37  23R:40  24R:37  25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln38         1R:37            2R:40         3R:40      4R:40  5R:40     6R:40    7R:37    8R:40   9R:40     10R:37      11R:40     12R:40  13R:40  14R:40  15R:40  16R:33  17R:37  18R:37  19R:37  20R:37  21R:33  22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln39         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12R:33  13R:37  14R:37  15R:37  16R:37  17R:37  18R:40  19R:40  20R:40  21R:40  22R:40  23R:40  24R:40  25R:40  26R:40  27R:37  28R:40  29R:40  30R:40  31R:40  32R:37
|aln40         1R:37            2R:40         3R:40      4R:40  5R:40     6R:40    7R:37    8R:40   9R:33     10R:40      11R:40     12R:37  13R:40  14R:40  15R:40  16R:37  17R:37  18R:37  19R:33  20R:37  21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln41         1R:40            2R:40         3R:40      4R:37  5R:40     6R:37    7R:40    8R:40   9R:40     10R:40      11R:40     12R:40  13R:40  14R:40  15R:40  16R:37  17R:40  18R:40  19R:40  20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln42         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21G:37  22R:40  23R:37  24R:37  25R:40  26R:27  27R:37  28R:40  29R:40  30R:40  31R:37  32R:27
|aln43         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21R:37  22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:37  30R:27  31R:33  32R:37
|aln44         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21S:4   22R:37  23R:40  24R:40  25R:37  26R:40  27R:40  28R:40  29R:37  30R:40  31R:40  32R:40
|aln45         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21S:7   22R:40  23R:40  24R:40  25R:40  26R:40  27R:40  28R:40  29R:40  30R:40  31R:40  32R:40
|aln46         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21S:7   22R:40  23R:40  24R:37  25R:37  26R:37  27R:15  28R:37  29R:40  30R:37  31R:33  32R:40
|aln47         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:37    8R:37   9R:37     10R:33      11R:33     12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln48         1R:40            2R:40         3R:40      4R:40  5R:40     6R:40    7R:37    8R:37   9R:37     10R:37      11R:33     12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln49         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24R:37  25R:37  26R:37  27R:37  28R:33  29R:37  30R:40  31R:37  32R:33
|aln50         1R:40            2R:40         3R:40      4R:37  5R:37     6R:37    7R:37    8R:33   9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln51         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25R:37  26R:37  27R:37  28R:37  29R:37  30R:40  31R:40  32R:40
|aln52         1R:40            2R:37         3R:33      4R:33  5R:40     6R:40    7R:33    8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln53         1R:37            2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32M:1
|aln54         1M:1             2M:1          3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11M:1      12M:1   13M:1   14M:1   15M:1   16M:1   17M:1   18M:1   19M:1   20M:1   21M:1   22M:1   23M:1   24M:1   25M:1   26M:1   27M:1   28M:1   29M:1   30M:1   31M:1   32R:27
|properties0   mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties1   mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties2   mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties3   mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties4   mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties5   mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties6   mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties7   mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties8   mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties9   mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties10  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties11  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties12  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties13  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties14  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties15  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties16  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties17  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties18  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties19  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties20  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties21  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties22  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties23  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties24  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties25  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties26  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties27  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties28  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties29  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties30  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties31  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties32  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties33  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties34  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties35  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties36  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties37  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties38  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties39  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties40  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties41  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties42  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties43  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties44  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties45  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties46  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties47  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties48  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties49  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties50  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1
|properties51  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties52  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties53  mapqual:60       strand:0      ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|properties54  mapqual:60       strand:1      ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1
|match0        1H:1             2H:0.2
|match1        1H:1             2H:0.2
|match2        1H:1             2H:0.2
|match3        1H:1             2H:0.2
|match4        1H:0.2           2H:1
|match5        1H:0.2           2H:1
|match6        1H:0.2           2H:1
|match7        1H:0.2           2H:1
|match8        1H:0.2           2H:1
|match9        1H:0.2           2H:1
|match10       1H:0.2           2H:1
|match11       1H:1             2H:0.2
|match12       1H:0.2           2H:1
|match13       1H:0.2           2H:1
|match14       1H:0.2           2H:1
|match15       1H:0.2           2H:1
|match16       1H:0.2           2H:1
|match17       1H:1             2H:0.2
|match18       1H:1             2H:0.2
|match19       1H:0.2           2H:1
|match20       1H:0.2           2H:1
|match21       1H:1             2H:0.2
|match22       1H:0.2           2H:1
|match23       1H:1             2H:0.2
|match24       1H:1             2H:0.2
|match25       1H:0.2           2H:1
|match26       1H:1             2H:0.2
|match27       1H:1             2H:0.2
|match28       1H:1             2H:0.2
|match29       1H:1             2H:0.2
|match30       1H:0.2           2H:1
|match31       1H:0.2           2H:1
|match32       1H:1             2H:0.2
|match33       1H:0.2           2H:1
|match34       1H:1             2H:0.2
|match35       1H:1             2H:0.2
|match36       1H:1             2H:0.2
|match37       1H:1             2H:0.2
|match38       1H:1             2H:0.2
|match39       1H:1             2H:0.2
|match40       1H:0.8           2H:0.2
|match41       1H:0.6           2H:0.2
|match42       1H:0             2H:0
|match43       1H:0.2           2H:0
|match44       1H:0             2H:0
|match45       1H:0             2H:0
|match46       1H:0             2H:0
|match47       1H:0             2H:0
|match48       1H:0             2H:0
|match49       1H:0             2H:0
|match50       1H:0             2H:0
|match51       1H:0             2H:0
|match52       1H:0             2H:0
|match53       1H:0             2H:0
|match54       1H:0             2H:0
|software      AC_1:1           QUAL:359.947

```

The first entry in the line defines the class of the example. By convention, we say the example is 1 if the haplotypes are correct, and -1 otherwise.

* `ref` : the reference
* `hap*`: the haplotypes described in the VCF record (including reference allele)
* `geno*`: the genotypes described in the VCF record (e.g. we have 2N of these for dipliods)
* `aln*` : the alignments overlapping the locus
* `properties*` : properties of the alignment taken from the input BAM file
* `match*` : the pairwise matching match of the alignment and each haplotype
* `software`: annotations from the VCF file, often software specific features

The reference sequence is given with the underlying bases in the `|ref` namespace.
    
The properties namespaces have these features:
    
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

