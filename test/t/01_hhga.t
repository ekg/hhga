#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 4

hhga -h 2>/dev/null
is $? 0 "hhga help runs"

is $(hhga -b minigiab/NA12878.chr22.tiny.bam -f minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -t -r q:10502-10562 | md5sum | cut -f 1 -d\ ) 8515e1483fbdd2121523c40310a8d324 "expected output produced for a test region"

is $(hhga -b minigiab/NA12878.chr22.tiny.bam -f minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -r q:10502-10562 -c 1 | md5sum | cut -f 1 -d\ ) 65464c66076ac7aa018822bad57844b5 "expected vw-format output produced for a test region"

is $(hhga -b minigiab/NA12878.chr22.tiny.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -t | grep ^hap | grep '\.----' | wc -l ) 1 "a normalized left-aligned indel is properly handled in the haplotypes"
