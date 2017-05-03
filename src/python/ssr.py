
from __future__ import division
import re
from Bio import SeqIO


def find_rep(seq_string, regex_dict, smaller_dict):
    # def sort_keys(key, rep_sub_dict):
    #     return (rep_sub_dict[key][0], key)

    for rep_length in regex_dict:
        rep_dict = {}
        for match in re.finditer(regex_dict[rep_length],
                                 seq_string):
            rep_seq, single_rep = match.groups()
            smaller_list = []
            for smaller_length in smaller_dict[rep_length]:
                regex = smaller_dict[rep_length][smaller_length]
                smaller_list.extend(re.findall(regex, single_rep))
            if len(smaller_list) == 0:
                rep_count = len(rep_seq) // len(single_rep)
                rep_dict[single_rep] = (rep_count, len(rep_seq), match.start())
        if rep_dict != {}:
            yield rep_length, rep_dict
        # sorted_rep = sorted(rep_dict,
        #                     key=lambda elem: sort_keys(elem, rep_dict),
        #                     reverse=True)
        # i = 0
        # while i < len(sorted_rep):
        #     for j in range(i + 1, len(sorted_rep)):
        #         this_rep = sorted_rep[i]
        #         other_rep = sorted_rep[j]
        #         rep_mut = ['%s%s' % (this_rep[k:], this_rep[:k])
        #                    for k in range(-1, -rep_length, -1)]
        #         if any(other_rep == rep for rep in rep_mut):
        #             del sorted_rep[j]
        #             rep_dict[]
        #     i += 1


def run(min_rep=3, min_rep_length=1, min_match_length=33, seq_length=66):
    regex_dict = {}
    smaller_dict = {}
    for rep_length in range(min_rep_length, seq_length // min_rep):
        print "(?=.{%i,}$)(([gatc]{%i})\\2{%i,})" % (min_match_length,
                                                     rep_length,
                                                     min_rep - 1)
        this_min_rep = max(min_rep, min_match_length // rep_length)
        regex_dict[rep_length] = re.compile("(([gatc]{%i})\\2{%i,})"
                                            % (rep_length,
                                               this_min_rep - 1),
                                            re.IGNORECASE)
        smaller_dict[rep_length] = {}
        for smaller_length in range(1, (rep_length // 2 + 1)):
            regex = re.compile("(([gatc]{%i})\\2{%i})"
                               % (smaller_length,
                                  rep_length // smaller_length - 1),
                               re.IGNORECASE)
            smaller_dict[rep_length][smaller_length] = regex
    bc_regex = re.compile(r'_([ACGT]+)')
    rep_dict = {}
    for seq in SeqIO.parse(('/home/cleemans/SURFdrive/TRIP/workspace/'
                            'mapping_table_rev.txt'), 'fasta'):
        this_seq_rep = find_rep(str(seq.seq), regex_dict, smaller_dict)
        bc = re.search(bc_regex, seq.id).group(1)
        for rep_length, new_dict in this_seq_rep:
            if rep_length in rep_dict:
                rep_bc_dict = rep_dict[rep_length]
                if bc in rep_bc_dict:
                    rep_single_dict = rep_bc_dict[bc]
                    for single_rep in new_dict:
                        if single_rep in rep_single_dict:
                            this_count = new_dict[single_rep]
                            count_dict = rep_single_dict[single_rep]
                            if this_count in count_dict:
                                count_dict[this_count] += 1
                            else:
                                count_dict[this_count] = 1
                        else:
                            rep_single_dict[single_rep] = {this_count: 1}
                else:
                    print {single_rep: {count: 1} for single_rep, count
                           in new_dict.items()}
                    rep_bc_dict[bc] = {single_rep: {count: 1}
                                       for single_rep, count
                                       in new_dict.items()}
            else:
                rep_dict[rep_length] = {bc: {single_rep: {count: 1}
                                             for single_rep, count
                                             in new_dict.items()}}
    return rep_dict
