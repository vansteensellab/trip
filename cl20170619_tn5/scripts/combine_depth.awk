#!/usr/bin/awk -f
{
    if (NR==1){
        chr=$1
        start=$2
        end=$3
        max_score=$4
        max_pos=1
        l=1
    } else {
        if ($1==chr && $2 == end){
            l += 1
            end = $3
            max_pos = $4>max_score?l:max_pos
            max_score = $4>max_score?$4:max_score
        } else {
            ori = max_pos > (l - max_pos)?"-":"+"
            print chr"\t"start"\t"end"\t"ori"\t"max_score
            chr=$1
            start=$2
            end=$3
            max_score=$4
            max_pos=1
            l=1
        }
    }
}

END {
    ori = max_pos > (length - max_pos)?"-":"+"
    print chr"\t"start"\t"end"\t"ori"\t"max_score
}
