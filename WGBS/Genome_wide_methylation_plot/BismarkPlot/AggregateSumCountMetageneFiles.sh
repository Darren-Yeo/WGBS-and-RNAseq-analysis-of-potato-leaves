
#!/bin/bash
##summarise by fragments (aggregate sum/aggregate counts)
path=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE

cd $path

for TE_files in Metagene_TE_*;do

echo ${TE_files}

awk '
BEGIN { FS=OFS="\t" }
NR==1 { print "fragment", "context", "meth_level"; next }

{
  key = $7 FS $5
  sum_sum[key] += $8
  sum_count[key] += $9
}

END {
  for (k in sum_sum) {
    meth = (sum_count[k] > 0) ? sum_sum[k] / sum_count[k] : 0
    print k, meth
  }
}
' ${TE_files} | (head -n 1 && tail -n +2 | sort -k1,1n) > Summary_${TE_files}

done

exit