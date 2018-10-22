#!/usr/bin/awk -f
BEGIN {
    OFS="\t"
    print "##gff-version 3"
}
{
    if ($1 ~ /^[^#]/){
        print $1, "tn5", "transgenic_insertion", $2, $3, ".", $4, ".", \
              "ID="name":"$1"_"$2"_"$5";Parent="name";Allele="$5 \
              ";CAST_mut="$6";129S1_mut="$7";depth_fwd="$8";depth_rev="$9 \
              ";avg_mapq="$10";max_mapq="$11
    }
}
