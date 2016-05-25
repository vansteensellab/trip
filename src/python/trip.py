#!/usr/bin/python
from __future__ import division
from collections import Counter
import re
from Bio.Data import IUPACData
from Bio import SeqIO
from itertools import izip
import Levenshtein
import warnings
import os
import pysam
import argparse
import itertools
from multiprocessing import Pool


class BKFilteredTree:
    def __init__(self, distfn, sorted_words, n, **kwargs):
        """
        Create a new BK-tree from the given distance function and
        words.

        Arguments:
        distfn: a binary function that returns the distance between
        two words.  Return value is a non-negative integer.  the
        distance function must be a metric space.

        words: an iterable.  produces values that can be passed to
        distfn

        """
        self.distfn = distfn
        self.kwargs = kwargs
        it = iter(sorted_words)
        root = it.next()
        self.tree = (root, {})
        self.n = n
        j = 0
        for i in it:
            self._add_word(self.tree, i)
            j += 1

    def _add_word(self, parent, word):
        pword, children = parent
        d = self.distfn(word, pword, **self.kwargs)
        iList = [i for i in range(d - self.n, d) if i in children]
        iList.extend(i for i in range(d + 1, d + self.n + 1) if i in children)
        j = 0
        found = False
        while j < len(iList) and not found:
            found = self._has_close_item(children.get(iList[j]), word)
            j += 1

        if d in children and not found:
            self._add_word(children[d], word)
        elif not found:
            children[d] = (word, {})

    def _has_close_item(self, parent, word):
        """
        Return all words in the tree that are within a distance of `n'
        from `word`.
        Arguments:

        parent: the parent tree to start from
        word: a word to query on
        n: a non-negative integer that specifies the allowed distance
        from the query word.
        """
        def rec(parent):
            pword, children = parent
            found = False
            d = self.distfn(word, pword, **self.kwargs)
            if d <= self.n:
                return True
            i = d - self.n
            while i <= d + self.n and not found:
                child = children.get(i)
                if child is not None:
                    found = rec(child)
                i += 1
            return found

        return rec(parent)

    def add_word(self, word):
        self._add_word(self.tree, word)

    def has_close_item(self, word):
        return self._has_close_item(self, self.tree, word)

    def query(self, word, n):
        """
        Return all words in the tree that are within a distance of `n'
        from `word`.
        Arguments:

        word: a word to query on
        n: a non-negative integer that specifies the allowed distance
        from the query word.

        Return value is a list of tuples (distance, word), sorted in
        ascending order of distance.
        """
        def rec(parent):
            pword, children = parent
            d = self.distfn(word, pword, **self.kwargs)
            results = []
            if d <= n:
                results.append((d, pword))
            for i in range(d - n, d + n + 1):
                child = children.get(i)
                if child is not None:
                    results.extend(rec(child))
            return results

        # sort by distance
        return sorted(rec(self.tree))

    def keys(self):
        def rec(parent):
            pword, children = parent
            yield pword
            for i in children.keys():
                child = children.get(i)
                if child is not None:
                    for word in rec(child):
                        yield word
        return rec(self.tree)


class BKTree:
    def __init__(self, distfn, words, **kwargs):
        """
        Create a new BK-tree from the given distance function and
        words.

        Arguments:
        distfn: a binary function that returns the distance between
        two words.  Return value is a non-negative integer.  the
        distance function must be a metric space.

        words: an iterable.  produces values that can be passed to
        distfn

        """
        self.distfn = distfn
        self.kwargs = kwargs
        it = iter(words)
        root = it.next()
        self.tree = (root, {})
        for i in it:
            self._add_word(self.tree, i)

    def add_word(self, word):
        self._add_word(self.tree, word)

    def _add_word(self, parent, word):
        pword, children = parent
        d = self.distfn(word, pword, **self.kwargs)
        if d in children:
            self._add_word(children[d], word)
        else:
            children[d] = (word, {})

    def query(self, word, n):
        """
        Return all words in the tree that are within a distance of `n'
        from `word`.
        Arguments:

        word: a word to query on
        n: a non-negative integer that specifies the allowed distance
        from the query word.

        Return value is a list of tuples (distance, word), sorted in
        ascending order of distance.
        """
        def rec(parent):
            pword, children = parent
            d = self.distfn(word, pword, **self.kwargs)
            results = []
            if d <= n:
                results.append((d, pword))
            for i in range(d - n, d + n + 1):
                child = children.get(i)
                if child is not None:
                    results.extend(rec(child))
            return results

        # sort by distance
        return sorted(rec(self.tree))

    def keys(self):
        def rec(parent):
            pword, children = parent
            yield pword
            for i in children.keys():
                child = children.get(i)
                if child is not None:
                    for word in rec(child):
                        yield word
        return rec(self.tree)


##############################################################################
############************  Function bc_extract_norm  ****************##########
##############################################################################
# DESCRIPTION:
# A subroutine to read the fastq file (typically the Single Read Illumina
# sequencing runs for normalization and expression counts). It reads the file
# and then extracts the barcode sequence based on the given arguments. The
# barcodes are stored in a hash as keys. The value for each key (barcode) is an
# anonymous array.
# If the reads are normalization, then the first element contains the
# normalization counts of that barcode. If the reads are expression than
# expression counts of that barcode are added as the second element of this
# anonymous array.
#
# INPUT:
# It takes 8 arguments
#   a) $fastq -> The name of the fastq file to be read (with full path)
#   b) $type -> Either "norm" (for normalization reads) or "exp" for expression
#                reads.
#   c) $ind_len -> The length of the index sequence. Zero if no index is part
#                   of the read. (It is important to note that nothing is done
#                   with index sequence per se. This is just to get the
#                   coordinates of pattern match right.
#   d) $bc_len -> The length of the barcode.
#   e) $pat1 -> The sequence of the first constant part.
#   f) $pat2 -> The sequence of the second constant part.
#   g) $BCs -> The reference to the empty hash (or hash with gDNA counts if
#               the type is cDNA) in which the barcodes and their counts are
#               stored.
#   h) $idx -> the position of the anonymous array in the hash where to add the
#                count. For $type eq "norm" it should be 0 and for others it
#                should be 1 ..
#
# DEPENDS:
# It depends on
#   Hdist (see below)
#   make_regex (see below)
#
# WORKING:
# It first determines the length of the $pat1 and the $pat2 (constant parts
# of the read).
# It generates the regular expression for matching using subroutine make_regex.
# One example of such a regex is.
#       ^.{24,28}ACAACTCGAG(.{15,17})TGATCCTGCA.+
# One advantage of this regexp is that it allows for some (but not too much)
# flexibility in exact position of the match of the first pattern.
# It starts with the second line in the fastq file which is the first sequence
# read and matches the regexp pattern. If the pattern is matched then the
# barcode is extracted using $1 variable (it store the recently matched
#  pattern).
# If the regexp pattern is not matched then, it finds the Hamming distance
# using the Hamming_dist function between the $pat1 and the corresponding
# part of the read. If this is below (length $pat1)/7 then same procedure
# is repeated for $pat2 and the correponding part of the read. Because the
# variance of one nucleotide in the barcode length is allowed (see above)
# so the $pat2 is matched at three possible places in the read. If the Hdist
# is below (length $pat2)/7 then the barcode sequence is extracted.
# In the end the barcode sequence is fed as keys to the hash BCs.
# The value of each key (barcode)  is an anonymous array with first element
# is normalization count and the second one (if  the $type is "exp")
# is expression count. During the the procedure the subroutine also takes
# record of total reads and the reads which could be matched properly.
#
# OUTPUT:
# It returns the total reads, the reads matched and a Counter object (behaves
# as a dictionary with all barcodes counted.
###############################################################################


def bc_extract_count(fastq, bc_len, constant1='', constant2='',
                     ind_len=0):
    bc_iterator = bc_extract(fastq, bc_len, constant1, constant2, ind_len)
    bc_count_dict = Counter(bc_iterator)
    return bc_count_dict, sum(bc_count_dict.values())


def bc_extract(fastq_name, bc_len, constant1='',
               constant2='', ind_len=0):
    len1 = len(constant1)
    len2 = len(constant2)
    regex = make_regex('exp', bc_len, constant1, constant2, ind_len)
    with open(fastq_name) as fastq:
        for line in fastq.readlines()[1::4]:
            line = line.strip()
            match = re.match(regex, line)
            bc_this = None
            if match is not None:
                bc_this = match.group(1)
            else:
                hdist_1 = hamming_dist(line[ind_len:(ind_len + len1)],
                                       constant1)
                if hdist_1 <= (len1 / 7):
                    start2 = (ind_len + len1 + bc_len)
                    end2 = (start2 + len2)
                    this_constant2 = line[(start2 - 1):(end2 + 1)]
                    this_len2 = end2 - start2
                    if hamming_dist(this_constant2[1:-1],
                                    constant2) < (this_len2 / 7):
                        bc_this = line[(ind_len + len1):start2]
                        # print 'first:'
                        # print line[(ind_len + len1):start2]
                        # print line
                    elif hamming_dist(this_constant2[:-2],
                                      constant2) < (this_len2 / 7):
                        bc_this = line[(ind_len + len1):(start2 - 1)]
                        # print 'second:'
                        # print line[(ind_len + len1):start2 - 1]
                        # print line
                    elif hamming_dist(this_constant2[2:],
                                      constant2) < (this_len2 / 7):
                        bc_this = line[(ind_len + len1):(start2 + 1)]
                        # print 'third:'
                        # print line[(ind_len + len1):start2 + 1]
                        # print line
            if bc_this is not None:
                yield(bc_this)


def mapping_fwd(fastq_in, fastq_out, constant1, constant2, bc_len, ind_len,
                bc_count_dict):
    len1 = len(constant1)
    len2 = len(constant2)
    total = 0
    hits = 0
    genuine = 0
    regex = make_regex('map', bc_len, constant1, constant2, ind_len)
    fin = open(fastq_in)
    fout = open(fastq_out, 'w')
    id_dict = {}
    for seq in SeqIO.parse(fin, 'fastq'):
        total += 1
        match = re.match(regex, str(seq.seq))
        bc_this = None
        if match is not None:
            bc_this = match.group(1)
            new_seq = seq[-len(match.group(2)):]
        else:
            hdist_1 = hamming_dist(seq.seq[ind_len:len1],
                                   constant1)
            constant_len = ind_len + len1 + len2 + bc_len
            if hdist_1 <= (len1 / 7):
                start2 = (ind_len + len1 + bc_len)
                if hamming_dist(seq.seq[start2:len1],
                                constant2) < (len2 / 7):
                    bc_this = str(seq.seq[(ind_len + len1):start2])
                    new_seq = seq[constant_len:]
                elif hamming_dist(seq.seq[(start2 - 1):len1],
                                  constant2) < (len2 / 7):
                    bc_this = str(seq.seq[(ind_len + len1):(start2 - 1)])
                    new_seq = seq[(constant_len - 1):]
                elif hamming_dist(seq.seq[(start2 + 1):len1],
                                  constant2) < (len2 / 7):
                    bc_this = str(seq.seq[(ind_len + len1):(start2 + 1)])
                    new_seq = seq[(constant_len + 1):]
        if bc_this is not None:
            hits += 1
            if bc_this in bc_count_dict:
                genuine += 1
            id_dict[seq.id] = bc_this
            new_id = '_'.join([seq.id, bc_this])
            new_seq.description = re.sub(seq.id, new_id,
                                         seq.description)
            new_seq.id = new_id
            SeqIO.write(new_seq, fout, 'fastq')
    fin.close()
    fout.close()
    return(total, hits, genuine, id_dict)


def mapping_rev(fastq_in, fastq_out, constant, id_dict):
    hits = 0
    length = len(constant)
    regex = make_regex('map', 0, constant)
    fin = open(fastq_in)
    fout = open(fastq_out, 'w')
    for seq in SeqIO.parse(fin, 'fastq'):
        if seq.id in id_dict:
            match = re.match(regex, str(seq.seq))
            if match is not None:
                new_seq = seq[-len(match.group(1)):]
                SeqIO.write(new_seq, fout, 'fastq')
                hits += 1
            elif (hamming_dist(str(seq.seq[:length]), constant) <
                  (length / 7)):
                new_seq = seq[length:]
                new_id = '_'.join((seq.id, id_dict[seq.id]))
                new_seq.description = re.sub(seq.id, new_id,
                                             seq.description)
                new_seq.id = new_id
                SeqIO.write(new_seq, fout, 'fastq')
                hits += 1
    return(hits)


def mapping(fastq_fwd, fastq_rev, fout_fwd, fout_rev, bc_len,
            constant_fwd1='', constant_fwd2='', constant_rev='',
            ind_len=0):
    len_fwd1 = len(constant_fwd1)
    len_fwd2 = len(constant_fwd2)
    len_rev = len(constant_rev)
    hits = 0
    missing_fwd = 0
    missing_rev = 0
    regex_fwd = make_regex('map', bc_len, constant_fwd1, constant_fwd2,
                           ind_len)
    regex_rev = make_regex('map', 0, constant_rev)
    fq_fwd = open(fastq_fwd)
    fq_rev = open(fastq_rev)
    fo_fwd = open(fout_fwd, 'w')
    fo_rev = open(fout_rev, 'w')
    it = izip(SeqIO.parse(fq_fwd, 'fastq'), SeqIO.parse(fq_rev, 'fastq'))
    for fwd, rev in it:
        if (fwd.id != rev.id):
            raise ValueError("forward and reverse id are not equal:\n%s v.s %s"
                             % (fwd.id, rev.id))
        match = re.match(regex_fwd, str(fwd.seq))
        bc_this = None
        if match is not None:
            bc_this = match.group(1)
            new_fwd = fwd[-len(match.group(2)):]
        else:
            hdist_1 = hamming_dist(fwd.seq[ind_len:len_fwd1], constant_fwd1)
            constant_len = ind_len + len_fwd1 + len_fwd2 + bc_len
            if hdist_1 <= (len_fwd1 / 7):
                start2 = (ind_len + len_fwd1 + bc_len)
                if hamming_dist(fwd.seq[start2:len_fwd2],
                                constant_fwd2) <= (len_fwd2 / 7):
                    bc_this = str(fwd.seq[(ind_len + len_fwd1):start2])
                    new_fwd = fwd[constant_len:]
                elif hamming_dist(fwd.seq[(start2 - 1):len_fwd2],
                                  constant_fwd2) <= (len_fwd2 / 7):
                    bc_this = str(fwd.seq[(ind_len + len_fwd1):(start2 - 1)])
                    new_fwd = fwd[(constant_len - 1):]
                elif hamming_dist(fwd.seq[(start2 + 1):len_fwd2],
                                  constant_fwd2) <= (len_fwd2 / 7):
                    bc_this = str(fwd.seq[(ind_len + len_fwd2):(start2 + 1)])
                    new_fwd = fwd[(constant_len + 1):]
        if bc_this is not None:
            match = re.match(regex_rev, str(rev.seq))
            if match is not None:
                new_rev = rev[-len(match.group(1)):]
            elif (hamming_dist(str(rev.seq[:len_rev]), constant_rev) <
                  (len_rev / 7)):
                new_rev = rev[len_rev:]
            else:
                new_rev = None
                missing_rev += 1
            if new_rev is not None:
                hits += 1
                new_id = '_'.join([fwd.id, bc_this])
                new_fwd.description = \
                    new_rev.description = re.sub(new_fwd.id, new_id,
                                                 new_fwd.description)
                new_fwd.id = new_rev.id = new_id
                SeqIO.write(new_fwd, fo_fwd, 'fastq')
                SeqIO.write(new_rev, fo_rev, 'fastq')
        else:
            missing_fwd += 1
            match = re.match(regex_rev, str(rev.seq))
            if (match is None and hamming_dist(str(rev.seq[:len_rev]),
                                               constant_rev) > (len_rev / 7)):
                missing_rev += 1
    fq_fwd.close()
    fq_rev.close()
    fo_fwd.close()
    fo_rev.close()
    return (hits, missing_fwd, missing_rev)


def unambiguous_ratio(string):
    return sum(string.count(char) / len(string)
               for char in IUPACData.unambiguous_dna_letters)


def hamming_dist(string1, string2):
    min_len = min(len(string1), len(string2))
    # base_dict = IUPACData.ambiguous_dna_values
    # hdist = sum(not any(i in base_dict[ch2] for i in base_dict[ch1])
    #             for ch1, ch2 in
    #             zip(string1[:min_len], string2[:min_len]))
    hdist = sum(ch1 != ch2 for ch1, ch2 in zip(string1[:min_len],
                string2[:min_len]))
    return(hdist)


def make_regex(kind, bc_len, constant1='', constant2='', ind_len=0):
    try:
        if (constant1 == '' and constant2 == ''):
            raise ValueError('missing constant sequence')
    except ValueError:
        print 'needs at least one constant sequence before or after barcode'

    len1 = len(constant1)
    len2 = len(constant2)
    if len1 > 10:
        ind_len += len1 - 10
    constant1 = constant1[-10:]
    if len2 > 10:
        constant2 = constant2[:10]
    min1 = max(0, ind_len - 2)
    max1 = ind_len + 2

    if kind == 'exp':
        end = '.+'
    elif kind == 'map' and (len2 - 10) > 0:
        end = '.{%i}(.+)' % len2 - 10
    elif kind == 'map' and (len2 - 10) <= 0:
        end = '(.+)'

    if bc_len != 0:
        pattern = r'^.{%i,%i}%s(.{%i,%i})%s%s'
        regex = pattern % (min1, max1, constant1, bc_len - 1, bc_len + 1,
                           constant2, end)
    else:
        pattern = r'^.{%i,%i}%s%s%s'
        regex = pattern % (min1, max1, constant1, constant2, end)
    print regex
    return re.compile(regex)


def get_mapping_string(sam_line):
    if sam_line.is_reverse:
        pos = sam_line.reference_end
    else:
        pos = sam_line.reference_start
    mapping = ':'.join([sam_line.reference_name, sam_line.flag, pos])
    return mapping


def config_parse(file_name, map_style):
    def check_integer(integer, default):
        if integer not in option_dict:
            warnings.warn('%s is not set, defaults to %i' % (integer, default))
            option_dict[integer] = default
        elif not option_dict[integer].isdigit():
            warnings.warn('The %s given is not an integer, defaults to %i'
                          % (integer, default))
            option_dict[integer] = default
        else:
            option_dict[integer] = int(option_dict[integer])

    def check_pat(pat):
        if pat not in option_dict:
            raise NameError("Could not find any specification of %s\n" % pat)
        elif re.search('[^ACGT]', option_dict[pat]) is not None:
            raise ValueError('The %s contains non-DNA characters\n' % pat)

    option_dict = {}
    with open(file_name) as f:
        for line in f.readlines():
            line = re.sub('#.*', '', line)
            option_list = line.split('=')
            option_list = [item.strip() for item in option_list]
            if len(option_list) == 2:
                option_dict[option_list[0]] = option_list[1]

    for pat in ['pat1', 'pat2']:
        check_pat(pat)

    for integer, default in [('index_length', 0),
                             ('barcode_length', 16),
                             ('lev_dist', 2),
                             ('min_counts', 5),
                             ('cores', 1)]:
        check_integer(integer, default)

    if map_style != 'n':
        for pat in ['map_pat1', 'map_pat2']:
            check_pat(pat)

        for integer, default in [('max_dist_for', 500),
                                 ('max_dist_rev', 20)]:
            check_integer(integer, default)

        if 'bowtie_base' not in option_dict:
            raise NameError("Could not find any specification of bowtie_base")
        elif not os.path.exists('%s.1.bt2' % option_dict['bowtie_base']):
            raise OSError('The bowtie_base for alignment %s does '
                          'not exist\n' % option_dict['bowtie_base'])
    if map_style == 'b' or map_style == 'r':
        check_pat('map_pat_rev')
    return option_dict


def parse_sam(sam_file, max_soft_clip=5, min_first_match=10,
              remap_soft_clip=17):
    remap_list = []
    map_dict = {}
    no_map_dict = {}
    for line in pysam.view("-B", 'test.bam'):
        bc_this = line.query_name.split('_')[1]
        if line.is_reverse:
            start_pos = line.reference_end
        else:
            start_pos = line.reference_start
        mapping = ':'.join((line.reference_name, line.flag, start_pos))
        add_mapping = False
        if line.cigarstring == '*':
            add_mapping = True
        else:
            first_cigar = line.cigar_tuple[0]
            if first_cigar[0] == 4:
                if first_cigar[1] <= max_soft_clip:
                    add_mapping = True
                if first_cigar[1] >= remap_soft_clip:
                    quality_dict = {'phred_quality': line.query_qualities}
                    seq = SeqIO.SeqRecord(line.query_sequence,
                                          line.query_name,
                                          discription=line.query_name,
                                          letter_annotations=quality_dict)
                    remap_list.append(seq)
            if first_cigar[0] == 0 and first_cigar[1] >= 10:
                add_mapping = True
        if add_mapping:
            if bc_this in map_dict:
                if mapping in map_dict[bc_this]:
                    map_dict[bc_this][mapping][0] += 1
                    map_dict[bc_this][mapping][1] += line.mapping_quality
                else:
                    map_dict[bc_this][mapping] = [1, line.mapping_quality]
            else:
                map_dict[bc_this] = {mapping: [1, line.mapping_quality]}
        else:
            if bc_this in no_map_dict:
                if mapping in no_map_dict[bc_this]:
                    no_map_dict[bc_this][mapping][0] += 1
                    no_map_dict[bc_this][mapping][1] += line.mapping_quality
                else:
                    no_map_dict[bc_this][mapping] = [1, line.mapping_quality]
            else:
                no_map_dict[bc_this] = {mapping: [1, line.mapping_quality]}
    return(map_dict, remap_list, no_map_dict)


def get_arg_options():
    cmd_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)
    cmd_parser.add_argument('-n', '--normFile', dest='normalization_file',
                            help=('file containing normalization reads. '
                                  'Fastq (unzipped)'))
    cmd_parser.add_argument('-e', '--expFile', dest='expression_file',
                            help=('file containing expression reads. '
                                  'Fastq (unzipped)'))
    cmd_parser.add_argument('-f', '--mapFor', dest='mapping_forward',
                            help=('file containing forward mapping reads. '
                                  'Fastq (unzipped)'))
    cmd_parser.add_argument('-r', '--mapRev', dest='mapping_reverse',
                            help=('file containing reverse mapping reads. '
                                  'Fastq (unzipped)'))
    cmd_parser.add_argument('-c', '--config', dest='configuration_file',
                            help=('file containing extra arguments '
                                  '(for an example see below)'))
    cmd_parser.add_argument('-o', '--out_dir', dest='output_directory',
                            help=('path to the directory where output should '
                                  'be saved'))
    cmd_parser.add_argument('-m', '--map', dest='map_style',
                            help=('[n OR b OR r OR f]\n'
                                  'should mapping be done from both '
                                  'forward and reverse (b) '
                                  'or from forward reads only (f) '
                                  'or from reverse reads only (r).'
                                  'If (f) only --mapFor needs to be specified '
                                  'otherwise (b or r) both  --mapFor and '
                                  '--mapRev needs to be specified.'
                                  'The default value is undefined which means '
                                  'no mapping.'))
    cmd_parser.add_argument('-b', '--bcFile', dest='barcode_file',
                            help=('optional file containing barcodes. '
                                  'If this file is not provided, a list of '
                                  'barcodes is generated based on the '
                                  'normalization file(s)'))
    cmd_parser.add_argument('-v', '--verbose', action='store_true',
                            dest='verbose',
                            help=('Prints out running commentary during the '
                                  'execution of the script'))
    cmd_parser.add_argument('-d', '--debug', action='store_true',
                            dest='debug',
                            help=('Prints out the progress into a file '
                                  'debug_file.txt in the output directory'))
    cmd_parser.description = ('trippy version 0.0.1\n'
                              'Adapted from the trip tool from TRIP Analysis '
                              'Software Kit (TASK), by Waseem Akhtar, '
                              'Johann de Jong, Jelle ten Hoeve and Ludo Pagie '
                              '(w.akhtar\@nki.nl), 11sep2013'
                              'This program reads FASTQ files from TRIP '
                              'experiments, extracts out the barcodes, '
                              'filters out mutant barcodes and makes '
                              'a table of barcode counts in normalization and '
                              'expression reads. In addition it extracts the '
                              'genomic neighborhood sequences from the '
                              'mapping reads and finds the locations '
                              '(and their frequency) '
                              'associated with these barcodes.')
    cmd_parser.epilog = ('\nExample of a Configuration File:\n\n'
                         'index_length   = 10    # length of the index '
                         '(0 = if no index used)\n'
                         'barcode_length = 16    # length of the barcode\n'
                         'pat1           = GTCACAAGGGCCGGCCACAACTCGAG\n'
                         '                       # first constant part of '
                         'normalization and expression reads\n'
                         'pat2           = '
                         'TGATCCTGCAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTT\n'
                         '                       # second constant part of '
                         'normalization and expression reads\n'
                         'lev_dist       = 2     # Levenshtein distance '
                         'threshold to find barcode mutants\n'
                         'map_pat1       = GTCACAAGGGCCGGCCACAACTCGAG\n'
                         '                       # first constant part of'
                         'forward mapping read\n'
                         'map_pat2       = TGATC # '
                         'second constant part of forward mapping read\n'
                         'map_pat_rev    = '
                         'GTACGTCACAATATGATTATCTTTCTAGGGTTAA\n'
                         '                       # constant part of reverse '
                         'mapping read\n'
                         'cores          = 2     # the number of processors '
                         'to be used for Bowtie2 alignments\n'
                         'bowtie_base    = '
                         '/Users/wa/CoolShit/bowtie2-2.1.0/mm9\n'
                         '                       # the index of the genome to'
                         'align against\n'
                         'max_dist_for   = 500   # max distance to cluster '
                         'forward mapping positions\n'
                         'max_dist_rev   = 20    # max distance to cluster '
                         'reverse mapping positions\n'
                         'min_counts     = 5     # minimum number of counts '
                         'for a genuine barcode\n')
    return cmd_parser.parse_args()


def run_barcode_counter(file_list, bc_len, pat1, pat2, index_len, cores,
                        verbose):
    if cores > 1:
        pool = Pool(processes=cores)
        it = itertools.izip(file_list,
                            itertools.repeat(bc_len),
                            itertools.repeat(pat1),
                            itertools.repeat(pat2),
                            itertools.repeat(index_len),
                            itertools.repeat(verbose))
        out = pool.map(barcode_counter_star, it)
    else:
        out = []
        for i in range(0, len(file_list)):
            out.append(barcode_counter(file_list[i], bc_len, pat1, pat2,
                                       index_len, verbose))
    return out


def barcode_counter_star(variables):
    return barcode_counter(*variables)


def barcode_counter(file_name, bc_len, pat1, pat2, index_len, verbose=False):
    if verbose:
        print ('Counting barcodes in the file '
               '%s ....' % file_name)
    bc_dict, hits = bc_extract_count(file_name, bc_len, pat1, pat2, index_len)
    file_base = os.path.basename(file_name)
    if verbose:
        unique_length = len(bc_dict)
        print "Total number of unique barcodes in %s:\t%i" % (file_base,
                                                              unique_length)
    return (bc_dict, hits)


def parse_bc_file(barcode_file):
    with open(barcode_file) as bc_file:
        for line in bc_file.readlines():
            yield line.strip().split(' ')


def build_tree(bc_dict, min_counts):
    def sort_keys(key, bc_dict):
        return (unambiguous_ratio(key), sum(bc_dict[key].values()))
    sorted_count = sorted(bc_dict.keys(),
                          key=lambda elem: sort_keys(elem, bc_dict),
                          reverse=True)
    bc_tree = BKTree(Levenshtein.distance, (sorted_count[0]))
    bc_index = 1
    bc_dict[sorted_count[0]]['index'] = bc_index
    part_match_dict = {}
    exact_match_dict = {}
    for bc_this in sorted_count[1:]:
        near_list = bc_tree.query(bc_this, config_dict['lev_dist'])
        # if there is no barcode closeby in in Levenshtein distance already
        # in the tree to point to a mutation, assume it's genuine and
        # add it if it's above the threshold
        if (len(near_list) == 0 and (sum(bc_dict[bc_this].values() >=
                                         min_counts))):
            bc_tree.add_word(bc_this)
            this_count = bc_dict[bc_this]
            for key in this_count:
                if key in exact_match_dict:
                    exact_match_dict[key] += this_count[key]
                else:
                    exact_match_dict[key] = this_count[key]
            bc_index += 1
            this_count['index'] = bc_index
        # else if there is one barcode that is closeby in Levenshtein
        # distance, or there is a single barcode that is closest
        # assume the barcode was mutated from this 'genuine' barcode.
        # Add the count to the genuine barcode and add it to the sum
        # of partly matched barcodes
        elif len(near_list) == 1:
                # or (len(near_list) > 1 and
                #     near_list[0][0] < near_list[1][0]):
            this_count = bc_dict.pop(bc_this)
            other_bc = near_list[0][1]
            if other_bc in bc_dict:
                other_count = bc_dict[other_bc]
            else:
                bc_index += 1
                bc_dict[other_bc] = other_count = {'index': bc_index}
            for key in this_count:
                if key in part_match_dict:
                    part_match_dict[key] += this_count[key]
                else:
                    part_match_dict[key] = this_count[key]
                if key in other_count:
                    other_count[key] += this_count[key]
                else:
                    other_count[key] = this_count[key]
        else:
            del bc_dict[bc_this]
    return (bc_tree, bc_dict, part_match_dict, exact_match_dict)


def filter_barcodes_star(kwargs):
    return filter_barcodes(*kwargs)


def filter_barcodes(bc_tree, bc_count_dict, bc_set=None):
    exact_match = 0
    part_match = 0
    if bc_set is None:
        bc_set = set(bc_tree.keys())
    for bc_this in bc_count_dict.keys():
        # if there is a bc_dict already and barcode is in it, just count as
        # exact match
        if bc_this in bc_set:
            exact_match += bc_count_dict[bc_this]
        else:
            near_list = bc_tree.query(bc_this, config_dict['lev_dist'])
            if len(near_list) == 1:
                    # or (len(near_list) > 1 and
                    #     near_list[0][0] < near_list[1][0]):
                other_bc = near_list[0][1]
                part_match += 1
                if other_bc in bc_count_dict:
                    bc_count_dict[other_bc] += bc_count_dict.pop(bc_this)
                else:
                    bc_count_dict[other_bc] = bc_count_dict.pop(bc_this)
            else:
                del bc_count_dict[bc_this]
    return bc_count_dict, exact_match, part_match


def run_barcode_filter(bc_tree, bc_count_list, cores):
    bc_set = set(bc_tree.keys())
    if cores > 1:
        pool = Pool(processes=cores)
        it = itertools.izip(itertools.repeat(bc_tree),
                            bc_count_list,
                            itertools.repeat(bc_set))
        out = pool.map(filter_barcodes_star, it)
    else:
        out = []
        for i in range(0, len(file_list)):
            out.append(filter_barcodes(bc_tree, bc_count_list, bc_set))
    return out


if __name__ == '__main__':
    options = get_arg_options()
    config_file = options.configuration_file
    norm_file_list = options.normalization_file.split(',')
    exp_file_list = options.expression_file.split(',')
    out_dir = options.output_directory
    bc_file = options.barcode_file
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    map_style = options.map_style
    if map_style not in ['n', 'b', 'r', 'f']:
        raise ValueError('map_style is not recognized as a valid '
                         'mapping style, options are [n OR b OR r OR f]')
    elif map_style != 'n':
        map_file_for = options.mapping_forward
        if map_style != 'f':
            map_file_rev = options.mapping_reverse
    verbose = options.verbose
    print verbose
    if verbose:
        print "Parsing the configuration file ....."
    config_dict = config_parse(config_file, map_style)

    if options.debug:
        debug_file = open('/'.join((out_dir, 'debug_file.txt')), 'w')
        debug_file.write('Your script has started working!!!!!\n\n')
        debug_file.write('Here are the details of command line input:\n')
        debug_file.write('normFile(s):\t%s\n' % norm_file_list)
        debug_file.write('expFile(s):\t%s\n' % exp_file_list)
        debug_file.write('config:\t%s\n' % config_file)
        debug_file.write('outdir:\t%s\n' % out_dir)
        debug_file.write('map:\t\t%s\n' % map_style)
        if map_style != 'n':
            debug_file.write('mapFor:\t%s\n' % map_file_for)
            if map_style != 'f':
                debug_file.write('mapRev:\t%s\n' % map_file_rev)
        if verbose:
            debug_file.write('The verbose option is on')
        debug_file.write('\n\n')
        debug_file.write('Here is what I got from your config file:\n')
        debug_file.write('index_length:\t%i\n' % config_dict['index_length'])
        debug_file.write('barcode_length:\t%i\n' %
                         config_dict['barcode_length'])
        debug_file.write('pat1:\t%s\n' % config_dict['pat1'])
        debug_file.write('pat2:\t%s\n' % config_dict['pat2'])
        debug_file.write('lev_dist:\t%i\n' % config_dict['lev_dist'])
        if map_style != 'n':
            debug_file.write('map_pat1:\t%s\n' % config_dict['map_pat1'])
            debug_file.write('map_pat2:\t%s\n' % config_dict['map_pat2'])
            if map_style != 'f':
                debug_file.write('map_pat_rev:\t%s\n'
                                 % config_dict['map_pat_rev'])
            debug_file.write('cores:\t%i\n' % config_dict['cores'])
            debug_file.write('bowtie_base:\t%s\n' % config_dict['bowtie_base'])
            debug_file.write('max_dist_for:\t%i\n' %
                             config_dict['max_dist_for'])
            if map_style != 'f':
                debug_file.write('max_dist_rev:\t%i\n'
                                 % config_dict['max_dist_rev'])
        debug_file.write('min_counts:\t%i' % config_dict['min_counts'])

    file_list = []
    file_list.extend(norm_file_list)
    file_list.extend(exp_file_list)
    count_out = run_barcode_counter(file_list,
                                    config_dict['barcode_length'],
                                    config_dict['pat1'],
                                    config_dict['pat2'],
                                    config_dict['index_length'],
                                    config_dict['cores'],
                                    verbose)

    total_match_dict = {}
    exact_match_dict = {}
    part_match_dict = {}
    if bc_file is not None:
        bc_dict = {key: {'index': value}
                   for value, key in parse_bc_file(bc_file)}
        bc_tree = BKTree(Levenshtein.distance, bc_dict.keys())
        bc_count_list = [count_dict for count_dict, hits in count_out]
        name_list = [os.path.basename(file_name) for file_name in file_list]
        total_match_dict = {name: hits for name, (cd, hits) in zip(name_list,
                                                                   count_out)}
        exact_match_dict = {}
        part_match_dict = {}
        filter_out = run_barcode_filter(bc_tree, bc_count_list,
                                        config_dict['cores'])
        for name, (bc_count_dict, exact_match, part_match) in zip(name_list,
                                                                  filter_out):
            for bc_this in bc_count_dict:
                bc_dict[bc_this][name] = bc_count_dict[bc_this]
            exact_match_dict[name] = exact_match
            part_match_dict[name] = part_match
    else:
        bc_dict = {}
        for i in range(0, len(norm_file_list)):
            norm_base = os.path.basename(norm_file_list[i])
            bc_norm_dict, hits = count_out[i]
            total_match_dict[norm_base] = hits
            exact_match_dict[norm_base] = 0
            part_match_dict[norm_base] = 0
            for bc in bc_norm_dict:
                if bc in bc_dict:
                    bc_dict[bc][norm_base] = bc_norm_dict[bc]
                else:
                    bc_dict[bc] = {norm_base: bc_norm_dict[bc]}
        tree_out = build_tree(bc_dict, config_dict['min_counts'])
        bc_tree, bc_dict, part_match_dict, exact_match_dict = tree_out
        bc_count_list = [count_dict for count_dict, h in
                         count_out[len(norm_file_list):]]
        name_list = [os.path.basename(file_name)
                     for file_name in norm_file_list]
        filter_out = run_barcode_filter(bc_tree, bc_count_list,
                                        config_dict['cores'])
        for name, (bc_count_dict, exact_match, part_match) in zip(name_list,
                                                                  filter_out):
            for bc_this in bc_count_dict:
                bc_dict[bc_this][name] = bc_count_dict[bc_this]
            exact_match_dict[name] = exact_match
            part_match_dict[name] = part_match

    print total_match_dict
    print exact_match_dict
    print part_match_dict
    with open('/'.join((out_dir, 'bc_count.txt')), 'w') as f_out:
        first_line_list = ['bc']
        first_line_list.extend([os.path.basename(f) for f in file_list])
        f_out.write('\t'.join(first_line_list))
        f_out.write('\n')
        for bc in bc_dict:
            line_list = [bc]
            this_count_dict = bc_dict[bc]
            for file_name in first_line_list[1:]:
                if file_name in this_count_dict:
                    line_list.append(str(this_count_dict[file_name]))
                else:
                    line_list.append('0')
            f_out.write('\t'.join(line_list))
            f_out.write('\n')

    # ind_len = 10
    # bc_len = 16
    # min_count = 5
    # constant1 = 'GTCACAAGGGCCGGCCACAACTCGAG'
    # constant2 = 'TGATCCTGCAGTG'
    # map_pat1 = 'GTCACAAGGGCCGGCCACAACTCGAG'
    # map_pat2 = 'TGATC'
    # map_pat_rev = 'GTACGTCACAATATGATTATCTTTCTAGGGTT'
    # trip_folder = '/home/cleemans/SURFdrive/TRIP'
    # laura_folder = '/'.join(('/home/NFS/users/l.brueckner',
    #                          'TTRIP_K562',
    #                          'lb20160331_fastqs_TRIP_KRAB'))
    # fastq_base = '3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_13.fq'
    # fastq_name = '/'.join([laura_folder, fastq_base])
    # fastq_name = '/'.join([trip_folder, 'raw_data', fastq_base])
    # fwd_base = '3354_1_iPCR_laura_eva_altndx_R1_001_smplIdx_09.fastq'
    # rev_base = '3354_1_iPCR_laura_eva_altndx_R2_001_smplIdx_09.fastq'
    # fwd_fastq = '/'.join((trip_folder, 'raw_data', fwd_base))
    # rev_fastq = '/'.join((trip_folder, 'raw_data', rev_base))
    # bc_count_dict = bc_extract_exp1(fastq_name, bc_len, constant1, constant2,
    #                                 ind_len)
    # for_map_tbl = '/'.join([trip_folder, 'workspace', 'mapping_tabe_for.fq'])
    # rev_map_tbl = '/'.join([trip_folder, 'workspace',
    #                         'mapping_table_rev.fq.gz'])
    # mapping_fwd(fwd_fastq, for_map_tbl, map_pat1, map_pat2, bc_len, ind_len,
    #             bc_count_dict)

    # ratioDict = {}
    # for dna in countDict.keys():
    #     ratio = trip.unambiguous_ratio(dna)
    #     if ratio in ratioDict:
    #         ratioDict[ratio].append(dna)
    #     else:
    #         ratioDict[ratio] = [dna]
    if options.debug:
        debug_file.close()
