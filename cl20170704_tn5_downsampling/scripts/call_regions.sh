#!/usr/bin/env bash

bedtools closest -d -a $1 -b $2 | \
awk -v max_gap=$3 '{
    if ($11 <= max_gap && $11 != -1){
        print $0
    }
}'
