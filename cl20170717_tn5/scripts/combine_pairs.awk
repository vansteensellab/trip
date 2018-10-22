#!/usr/bin/awk -f
function abs(value)
{
  return (value<0?-value:value);
}


{
    if ($1 != "." && $4 != "."){
        if ($1 == $4){
            if (abs($3-$4) < max_dist || abs($5-$2) < max_dist){
                if ($9 != $10){
                    start=$2<$5?$2:$5
                    end=$3>$6?$3:$6
                    print $1"\t"start"\t"end"\t"$9
                    correct += 1
                } else {
                    same_ori += 1
                }
            } else {
                too_far += 1
            }
        } else {
            mapped_dif += 1
        }
    } else {
        if ($1 == $4){
            both_unmapped += 1
        } else {
            one_unmapped += 1
        }
    }
}

END {
    total=correct + same_ori + mapped_dif + both_unmapped + one_unmapped + too_far
    unmapped=both_unmapped + one_unmapped
    print "total number of read-pairs:\t         "total > "/dev/stderr"
    print "\tcorrectly mapped:                "correct > "/dev/stderr"
    print "\tpairs in same orientation:       "same_ori > "/dev/stderr"
    print "\tpairs are too far apart:         "too_far > "/dev/stderr"
    print "\tmapped to different chromosomes: "mapped_dif > "/dev/stderr"
    print "\treads not mapped:                "unmapped > "/dev/stderr"
    print "\t\tone read not mapped:     "one_unmapped > "/dev/stderr"
    print "\t\tboth reads not mapped:   "both_unmapped > "/dev/stderr"
}
