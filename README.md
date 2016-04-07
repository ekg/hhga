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

As an example, if this is a tview-like representation of a set of alignments at a putative SNP site, augmented with candidate haplotypes. Generated via `hhga -b minigiab/9251-9252.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -w 20 -t` from the `test/` directory.

```txt
q_9251_GTTCT_G
reference   ATTGTGCCAAGTTCTTTCTT
hap                   .....
hap                   .----
SodqfXYPZI  ..                   60 chr22.bin8.cram:166:8383
SodqfXYPZI  .....                60 chr22.bin8.cram:166:8410
sOdqfXYPZI  ..................   60 chr22.bin8.cram:166:8436
sOdqFxYPZI  ...........----..... 60 chr22.bin8.cram:166:8455
SodqFxYPZI  ...........----..... 60 chr22.bin8.cram:166:8390
sOdqFxYPZI  .................... 60 chr22.bin8.cram:166:8457
SodqFxYPZI  .................... 60 chr22.bin8.cram:166:8436
SodqfXYPZI  ...........----..... 60 chr22.bin8.cram:166:8400
sOdqfXYPZI  ...........----..... 60 chr22.bin8.cram:166:8460
sOdqfXYPZI       ............... 60 chr22.bin8.cram:166:8475
SodqFxYPZI                   ... 60 chr22.bin8.cram:166:8395
```    

We use six symbols to encode the MSA, `{ A, T, G, C, N, U, M, R }`, where `N` is the degenerate symbol (it could represent any of A, T, G, or C), `U` is a gap symbol required to normalize the MSA into a matrix, and `M` is a symbol indicating if we don't have any information at the position in the MSA. In standard practice, we use `R` to indicate when a base is the same as the reference base. This reduces the feature complexity of the model and helps performance. We further extend this to one specific symbol for each position in the matrix. So the reference-matching base at position 23 is always `23R`.

Various features of the reads are represented in other namespaces.

The feature space transformation of this site would be given by `hhga -b minigiab/9251-9252.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -w 20 -c 1`. The output will come on a single line.

```txt
1 'q_9251_GTTCT_G |ref 1A:1 2T:1 3T:1 4G:1 5T:1 6G:1 7C:1 8C:1 9A:1 10A:1 11G:1 12T:1 13T:1 14C:1 15T:1 16T:1 17T:1 18C:1 19T:1 20T:1 |hap0 1M:1 2M:1 3M:1 4M:1 5M:1 6M:1 7M:1 8M:1 9M:1 10M:1 11R:1 12R:1 13R:1 14R:1 15R:1 16M:1 17M:1 18M:1 19M:1 20M:1 |hap1 1M:1 2M:1 3M:1 4M:1 5M:1 6M:1 7M:1 8M:1 9M:1 10M:1 11R:1 12U:1 13U:1 14U:1 15U:1 16M:1 17M:1 18M:1 19M:1 20M:1 |aln0 mapqual:0.999999 strand:1 ostrand:0 dup:0 qcfail:0 fmate:0 xmate:1 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9998 2R:0.999499 3M:1 4M:1 5M:1 6M:1 7M:1 8M:1 9M:1 10M:1 11M:1 12M:1 13M:1 14M:1 15M:1 16M:1 17M:1 18M:1 19M:1 20M:1 |aln1 mapqual:0.999999 strand:1 ostrand:0 dup:0 qcfail:0 fmate:0 xmate:1 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9998 2R:0.9998 3R:0.9998 4R:0.999499 5R:0.999499 6M:1 7M:1 8M:1 9M:1 10M:1 11M:1 12M:1 13M:1 14M:1 15M:1 16M:1 17M:1 18M:1 19M:1 20M:1 |aln2 mapqual:0.999999 strand:0 ostrand:1 dup:0 qcfail:0 fmate:0 xmate:1 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9999 2R:0.9999 3R:0.99369 4R:0.9998 5R:0.9998 6R:0.9998 7R:0.9998 8R:0.9998 9R:0.9999 10R:0.9998 11R:0.9998 12R:0.998005 13R:0.999499 14R:0.9999 15R:0.9999 16R:0.9998 17R:0.9999 18R:0.9998 19M:1 20M:1 |aln3 mapqual:0.999999 strand:0 ostrand:1 dup:0 qcfail:0 fmate:1 xmate:0 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9999 2R:0.9999 3R:0.9999 4R:0.9999 5R:0.9999 6R:0.9999 7R:0.9999 8R:0.9999 9R:0.9999 10R:0.999499 11R:0.9999 12U:0.9999 13U:0.999499 14U:0.9999 15U:0.9999 16R:0.9999 17R:0.9999 18R:0.9999 19R:0.9999 20R:0.9999 |aln4 mapqual:0.999999 strand:1 ostrand:0 dup:0 qcfail:0 fmate:1 xmate:0 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9999 2R:0.998005 3R:0.9998 4R:0.9999 5R:0.9999 6R:0.9999 7R:0.9999 8R:0.9999 9R:0.9998 10R:0.9999 11R:0.9999 12U:0.9998 13U:0.9999 14U:0.9999 15U:0.9999 16R:0.9999 17R:0.9999 18R:0.9999 19R:0.9999 20R:0.9999 |aln5 mapqual:0.999999 strand:0 ostrand:1 dup:0 qcfail:0 fmate:1 xmate:0 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9999 2R:0.9999 3R:0.999499 4R:0.999499 5R:0.9999 6R:0.9999 7R:0.9999 8R:0.9999 9R:0.968377 10R:0.999499 11R:0.998005 12R:0.999499 13R:0.968377 14R:0.9998 15R:0.9999 16R:0.9999 17R:0.999499 18R:0.9999 19R:0.99369 20R:0.9999 |aln6 mapqual:0.999999 strand:1 ostrand:0 dup:0 qcfail:0 fmate:1 xmate:0 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9998 2R:0.999499 3R:0.9999 4R:0.9998 5R:0.9999 6R:0.9999 7R:0.9999 8R:0.9999 9R:0.9999 10R:0.9999 11R:0.9999 12R:0.9999 13R:0.9999 14R:0.9998 15R:0.9999 16R:0.9999 17R:0.9999 18R:0.9999 19R:0.9999 20R:0.9999 |aln7 mapqual:0.999999 strand:1 ostrand:0 dup:0 qcfail:0 fmate:0 xmate:1 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9999 2R:0.9999 3R:0.9998 4R:0.9998 5R:0.9999 6R:0.9998 7R:0.9999 8R:0.9998 9R:0.9998 10R:0.9999 11R:0.9998 12U:0.9998 13U:0.9999 14U:0.9998 15U:0.9999 16R:0.9999 17R:0.9998 18R:0.9999 19R:0.9998 20R:0.9999 |aln8 mapqual:0.999999 strand:0 ostrand:1 dup:0 qcfail:0 fmate:0 xmate:1 ymap:1 paired:1 zprimary:1 iproper:1 1R:0.9999 2R:0.9999 3R:0.9999 4R:0.9999 5R:0.9999 6R:0.9999 7R:0.9999 8R:0.9999 9R:0.9999 10R:0.9999 11R:0.9999 12U:0.9999 13U:0.9999 14U:0.9999 15U:0.9999 16R:0.9999 17R:0.9999 18R:0.9999 19R:0.9999 20R:0.9999 |aln9 mapqual:0.999999 strand:0 ostrand:1 dup:0 qcfail:0 fmate:0 xmate:1 ymap:1 paired:1 zprimary:1 iproper:1 1M:1 2M:1 3M:1 4M:1 5M:1 6R:0.999499 7R:0.9998 8R:0.9998 9R:0.9998 10R:0.9998 11R:0.9998 12R:0.9999 13R:0.9999 14R:0.9999 15R:0.9999 16R:0.9999 17R:0.9999 18R:0.9999 19R:0.9999 20R:0.9999 |aln10 mapqual:0.999999 strand:1 ostrand:0 dup:0 qcfail:0 fmate:1 xmate:0 ymap:1 paired:1 zprimary:1 iproper:1 1M:1 2M:1 3M:1 4M:1 5M:1 6M:1 7M:1 8M:1 9M:1 10M:1 11M:1 12M:1 13M:1 14M:1 15M:1 16M:1 17M:1 18R:0.9998 19R:0.9998 20R:0.9998
```

For debugging we can pipe the hhga format output through `sed s/\|/\\n\|/g | column -t` to convert the spaces into newlines and columnarize the fields.

```txt
1       'q_9251_GTTCT_G
    |ref    1A:1              2T:1      3T:1       4G:1   5T:1      6G:1     7C:1     8C:1    9A:1      10A:1       11G:1      12T:1      13T:1        14C:1        15T:1        16T:1        17T:1        18C:1      19T:1      20T:1
    |hap0   1M:1              2M:1      3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11R:1      12R:1      13R:1        14R:1        15R:1        16M:1        17M:1        18M:1      19M:1      20M:1
    |hap1   1M:1              2M:1      3M:1       4M:1   5M:1      6M:1     7M:1     8M:1    9M:1      10M:1       11R:1      12U:1      13U:1        14U:1        15U:1        16M:1        17M:1        18M:1      19M:1      20M:1
    |aln0   mapqual:0.999999  strand:1  ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9998  2R:0.999499  3M:1         4M:1         5M:1         6M:1         7M:1       8M:1       9M:1         10M:1         11M:1         12M:1         13M:1         14M:1       15M:1       16M:1       17M:1         18M:1       19M:1        20M:1
    |aln1   mapqual:0.999999  strand:1  ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9998  2R:0.9998    3R:0.9998    4R:0.999499  5R:0.999499  6M:1         7M:1       8M:1       9M:1         10M:1         11M:1         12M:1         13M:1         14M:1       15M:1       16M:1       17M:1         18M:1       19M:1        20M:1
    |aln2   mapqual:0.999999  strand:0  ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9999  2R:0.9999    3R:0.99369   4R:0.9998    5R:0.9998    6R:0.9998    7R:0.9998  8R:0.9998  9R:0.9999    10R:0.9998    11R:0.9998    12R:0.998005  13R:0.999499  14R:0.9999  15R:0.9999  16R:0.9998  17R:0.9999    18R:0.9998  19M:1        20M:1
    |aln3   mapqual:0.999999  strand:0  ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9999  2R:0.9999    3R:0.9999    4R:0.9999    5R:0.9999    6R:0.9999    7R:0.9999  8R:0.9999  9R:0.9999    10R:0.999499  11R:0.9999    12U:0.9999    13U:0.999499  14U:0.9999  15U:0.9999  16R:0.9999  17R:0.9999    18R:0.9999  19R:0.9999   20R:0.9999
    |aln4   mapqual:0.999999  strand:1  ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9999  2R:0.998005  3R:0.9998    4R:0.9999    5R:0.9999    6R:0.9999    7R:0.9999  8R:0.9999  9R:0.9998    10R:0.9999    11R:0.9999    12U:0.9998    13U:0.9999    14U:0.9999  15U:0.9999  16R:0.9999  17R:0.9999    18R:0.9999  19R:0.9999   20R:0.9999
    |aln5   mapqual:0.999999  strand:0  ostrand:1  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9999  2R:0.9999    3R:0.999499  4R:0.999499  5R:0.9999    6R:0.9999    7R:0.9999  8R:0.9999  9R:0.968377  10R:0.999499  11R:0.998005  12R:0.999499  13R:0.968377  14R:0.9998  15R:0.9999  16R:0.9999  17R:0.999499  18R:0.9999  19R:0.99369  20R:0.9999
    |aln6   mapqual:0.999999  strand:1  ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9998  2R:0.999499  3R:0.9999    4R:0.9998    5R:0.9999    6R:0.9999    7R:0.9999  8R:0.9999  9R:0.9999    10R:0.9999    11R:0.9999    12R:0.9999    13R:0.9999    14R:0.9998  15R:0.9999  16R:0.9999  17R:0.9999    18R:0.9999  19R:0.9999   20R:0.9999
    |aln7   mapqual:0.999999  strand:1  ostrand:0  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9999  2R:0.9999    3R:0.9998    4R:0.9998    5R:0.9999    6R:0.9998    7R:0.9999  8R:0.9998  9R:0.9998    10R:0.9999    11R:0.9998    12U:0.9998    13U:0.9999    14U:0.9998  15U:0.9999  16R:0.9999  17R:0.9998    18R:0.9999  19R:0.9998   20R:0.9999
    |aln8   mapqual:0.999999  strand:0  ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1  1R:0.9999  2R:0.9999    3R:0.9999    4R:0.9999    5R:0.9999    6R:0.9999    7R:0.9999  8R:0.9999  9R:0.9999    10R:0.9999    11R:0.9999    12U:0.9999    13U:0.9999    14U:0.9999  15U:0.9999  16R:0.9999  17R:0.9999    18R:0.9999  19R:0.9999   20R:0.9999
    |aln9   mapqual:0.999999  strand:0  ostrand:1  dup:0  qcfail:0  fmate:0  xmate:1  ymap:1  paired:1  zprimary:1  iproper:1  1M:1       2M:1         3M:1         4M:1         5M:1         6R:0.999499  7R:0.9998  8R:0.9998  9R:0.9998    10R:0.9998    11R:0.9998    12R:0.9999    13R:0.9999    14R:0.9999  15R:0.9999  16R:0.9999  17R:0.9999    18R:0.9999  19R:0.9999   20R:0.9999
    |aln10  mapqual:0.999999  strand:1  ostrand:0  dup:0  qcfail:0  fmate:1  xmate:0  ymap:1  paired:1  zprimary:1  iproper:1  1M:1       2M:1         3M:1         4M:1         5M:1         6M:1         7M:1       8M:1       9M:1         10M:1         11M:1         12M:1         13M:1         14M:1       15M:1       16M:1       17M:1         18R:0.9998  19R:0.9998   20R:0.9998
    |software freebayesLength: 0.9 
    
```

The first entry in the line defines the class of the example. By convention, we say the example is 1 if the haplotypes are correct, and -1 otherwise.

* `ref` : the reference
* `hap0`: the first haplotype
* `hap1`: the second haplotype
* `aln*` : the alignments
* `software`: freebayes or vg specific features, lnl:1


Fold each of the following that correspond to each alignment:
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

