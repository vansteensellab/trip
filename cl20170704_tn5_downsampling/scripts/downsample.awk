#!/usr/bin/awk -f

BEGIN {
    count=0
}
{
    if (NR%4==1){
        count += 1
    }
    if (count < num){
        if (NR%4==2 || NR%4==0){
            print substr($0,0,len)
        } else {
            print $0
        }
    }
}
