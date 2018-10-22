#!/usr/bin/awk -f

{
    if ($1 ~ /^[^#]/){
        match($8, "QS=([0-9.]*)",arr)
        if (arr[1] < 0.5){
            if (NR==FNR){
                m1 += 1
            } else {
                m2 += 1
            }
        }
    }
}
END {
    print m1"\t"m2
}
