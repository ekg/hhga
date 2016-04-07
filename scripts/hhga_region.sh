#!/bin/bash

if [ $# -ne 4 ];
then
echo usage: $0 '[window_size] [region] [vcf_dir] [output_dir]'
echo generates hhga data for this region in the current directory
exit
fi

region=$1
window_size=$2
vcf_dir=$3
output_dir=$4
ref=~/graphs/hs37d5.fa
dir=~/graphs/NA12878
input_bam=$dir/XPrize_Illumina_WG.bam
variant_calls=$vcf_dir/$region.vcf.gz
giab_truth=$dir/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
giab_callable=$dir/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz

giab_region=$output_dir/$region.giab_truth.vcf.gz
normalized_variants=$output_dir/$region.norm.vcf.gz
callable_variants=$output_dir/$region.giab_callable.vcf.gz
uncallable_variants=$output_dir/$region.giab_uncallable.vcf.gz
passing_variants=$output_dir/$region.giab_pass.vcf.gz
failing_variants=$output_dir/$region.giab_fail.vcf.gz
hhga_truthiness=$output_dir/$region.truthiness.hhga.gz

echo getting giab region
tabix -h $giab_truth $region | bgziptabix $giab_region

echo making normalized variants
tabix -h $variant_calls $region | vcffilter -f 'AC > 0' | vcfallelicprimitives -kg | vt normalize -r $ref -q - | bgziptabix $normalized_variants

echo intersecting with callable
bedtools intersect -header -a $normalized_variants -b $giab_callable | bgziptabix $callable_variants
bedtools intersect -header -v -a $normalized_variants -b $giab_callable | bgziptabix $uncallable_variants

echo intersecting with giab
vcfintersect -r $ref -i $giab_region $callable_variants | bgziptabix $passing_variants
vcfintersect -v -r $ref -i $giab_region $callable_variants | bgziptabix $failing_variants

echo generating truthiness hhga

( hhga -f $ref -r $region -b $input_bam -v $passing_variants -w $window_size -c 1
  hhga -f $ref -r $region -b $input_bam -v $failing_variants -w $window_size -c -1
) | shuf | gzip >$hhga_truthiness

echo cleaning up

rm -f $giab_region*
rm -f $normalized_variants*
rm -f $callable_variants*
rm -f $uncallable_variants*
rm -f $passing_variants*
rm -f $failing_variants*

