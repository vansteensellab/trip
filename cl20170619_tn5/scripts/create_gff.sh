#!/usr/bin/env bash







for input in "$@":
    arr=($(echo $input))
    name=${arr[0]}
    insertions=${arr[1]}
    barcode=${arr[2]}
    awk -vname=$name '{
            if (NR==FNR){
                if ($3 != "*"){
                    arr[$3-$4] = $1
                }
            } else {
                print "$2\t"name"\tclone_insert\tstart\tend\t.\t+\t.\tID=insert0001;BARCODE=;ALLELE="
            }
        }'
