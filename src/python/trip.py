#!/home/NFS/users/c.leemans/python/bin/python2.7
from __future__ import division
from collections import Counter
import read_parser
import re
import gzip
from Bio.Data import IUPACData
from Bio import SeqIO
from Bio import Seq
import warnings
import os
import pysam
import argparse
import itertools
import subprocess
from multiprocessing import Pool
from datetime import datetime


##############################################################################
##########*******************  #Parsing reads#  ********************##########
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


def build_exp_structure(index_len, bc_len, pat1, pat2, structure_file):
    """ wrapper around create structure expression and normalization files.

    Typical structure of Normalization and expression reads:
       index            pat1               barcode           pat2
    NNNNNNNNNNGTCACAAGGGCCGGCCACAACTCGAGNNNNNNNNNNNNNNNNTGATCCTGCAGTGTCACC...

    Keyword arguments:
    index_len -- the length of the index at the beginning of the read
    bc_len -- the length of the barcode
    pat1 -- the first constant pattern (after index, before barcode)
    pat2 -- the pattern after the barcode
    sturcture_file -- output file to save the structure
    """
    exp_structure = structure(structure_file)
    if index_len > 0:
        exp_structure.add_segment('index', 'const', index_len)
    exp_structure.add_segment('pat1', 'const', pat1)
    exp_structure.add_segment('barcode', 'barcode', bc_len)
    exp_structure.add_segment('pat2', 'const', pat2, is_fixed=False)
    exp_structure.write()
    return(structure_file)


def build_mapping_structure(index_len, bc_len, map_pat1, map_pat2, map_pat_rev,
                            structure_file):
    """ wrapper around parse reads for forward and reverse mapping files.

    Typical structure of forward reads:
       index          map_pat1             barcode   map_pat2    gDNA
    NNNNNNNNNNGTCACAAGGGCCGGCCACAACTCGAGNNNNNNNNNNNNNNNNTGATCNNNNNNNNNNNN...

    Typical structure of reverse reads:
            map_pat_rev                             gDNA
    GTACGTCACAATATGATTATCTTTCTAGGGTTAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...

    Keyword arguments:
    index_len -- the length of the index at the beginning of the read
    bc_len -- the length of the barcode
    map_pat1 -- the pattern in between index and barcode on the forward read
    map_pat2 -- the pattern in between the barcode and the gDNA on the forward
                read
    map_pat_rev -- the pattern on the reverse read
    sturcture_file -- output file to save the structure
    """
    map_structure = structure(structure_file)
    print 'parsing mapping reads'
    if index_len > 0:
        map_structure.add_segment('index', 'const', index_len)

    # because it is not easy to locate map_pat2 in between two unknown sections
    # there is a trick in the read parser in which the barcode surrounded by
    # constant parts is chopped of first and subsequently dissected
    # into the barcode and the two constant parts.
    # This code corrects the error threshold to account for the unknown barcode
    # piece in the middle, but it does not care where the errors occur.
    # This means that the second pattern has all te errors, this could make for
    # situations where the construct passes the thresholds, but in discecting
    # the pieces, thresholds are not passed.
    # To fix that problem the first pattern can be split so a small part of it
    # is added to the barcode construct.
    if len(map_pat1) / len(map_pat2) > 1.5:
        map_pat1a = map_pat1[:-len(map_pat2)]
        const_bar = '%s{%i}%s' % (map_pat1[-len(map_pat2):], bc_len,
                                  map_pat2)
        map_structure.add_segment('map_pat1a', 'const', map_pat1a)
        map_structure.add_segment('const_bar', 'const_bar', const_bar)
    # hypothetically this could also work the other way around
    elif len(map_pat2) / len(map_pat1) > 1.5:
        const_bar = '%s{%i}%s' % (map_pat1, bc_len,
                                  map_pat2[:len(map_pat1)])
        map_pat2b = map_pat2[len(map_pat1):]
        map_structure.add_segment('const_bar', 'const_bar', const_bar)
        map_structure.add_segment('map_pat2b', 'const', map_pat2b)
    else:
        const_bar = '%s{%i}%s' % (map_pat1, bc_len, map_pat2)
        map_structure.add_segment('const_bar', 'const_bar', const_bar)
    # sometimes the iPCR product is really small and part of the read is
    # reading back into the constant part at the 3' end
    rev_comp = str(Seq.Seq(map_pat_rev).reverse_complement())

    # in the case of the reverse reads, the barcode which we don't know now
    # will be also in the 3' of really short iPCR products.
    # this barcode needs to be added later in the read parser.
    map_pat1_comp = str(Seq.Seq(map_pat1).reverse_complement())
    map_pat2_comp = str(Seq.Seq(map_pat2).reverse_complement())
    constant_bar_comp = '%s[BC]%s' % (map_pat2_comp, map_pat1_comp)
    map_structure.add_segment('rev_map_complement', 'const', rev_comp,
                              is_5=False, is_present=False, is_fixed=False)
    map_structure.add_segment('rev_map', 'const', map_pat_rev, is_second=True)
    map_structure.add_segment('fwd_map_complement', 'const_bar_comp',
                              constant_bar_comp, is_5=False, is_present=False,
                              is_second=True, is_fixed=False)
    map_structure.write()


class structure:
    def __init__(self, structure_file):
        self.structure_file = structure_file
        self.line_list = []

    def add_segment(self, part_id, part_type, dna, is_5=True,
                    is_present=True, is_second=False, is_fixed=True):
        if type(dna) == int:
            dna = str(dna)
        second = 'True' if is_second else 'False'
        pos = 'fixed' if is_fixed else 'var'
        present = 'present' if is_present else '-'
        if is_5:
            dna_5 = dna
            dna_3 = '-'
        else:
            dna_5 = '-'
            dna_3 = dna
        line = '\t'.join((part_id, dna_5, dna_3, part_type, present,
                          second, pos))
        self.line_list.append(line)

    def write(self):
        with open(self.structure_file, 'w') as f_out:
            f_out.write('\t'.join(("ID", "5'", "3'", "type", 'req',
                                   'second-read', 'pos')))
            f_out.write('\n')
            for line in self.line_list:
                f_out.write(line)
                f_out.write('\n')


def parse_reads(fastq_name, structure_file, out_dir):
    """ run the read parser and return barcode counts, statistics and the name
    of the file. Can be used with pair of forward and reverse read files, or
    with single fastq files.
    """
    if type(fastq_name) == str:
        base_name = get_base_noext(fastq_name)
        out_base = '/'.join((out_dir, base_name))
        fastq_name = [fastq_name, False]
    else:
        base_pair = [get_base_noext(f) for f in fastq_name]
        base_name = ''.join((ca for ca, cb in
                             zip(base_pair[0], base_pair[1])
                             if ca == cb))
        out_base = '/'.join((out_dir, base_name))
    read_parser.run(fastq_name, structure_file, out_base)
    bc_iter = parse_barcodes(''.join((out_base, '.barcode.txt.gz')))
    bc_count_dict = Counter(bc_vec[1] for bc_vec in bc_iter)
    stats = parse_parsing_stats(''.join((out_base, '.statistics.txt')))
    return(bc_count_dict, stats, base_name)


def parse_barcodes(barcode_file, file_type='gzip'):
    """ parse a barcode file produced by the read parser
    """
    if file_type == 'gzip':
        open_func = gzip.open
    else:
        open_func = open
    with open_func(barcode_file) as bc_file:
        bc_file.next()
        for line in bc_file.readlines():
            yield(line.strip().split('\t'))


def parse_parsing_stats(stat_file):
    with open(stat_file) as f:
        line_1 = f.readline()
        line_2 = f.readline()
    line_split_1 = line_1.strip().split('\t')
    line_split_2 = line_2.strip().split('\t')
    return({a: b for a, b in zip(line_split_1, line_split_2)})


def config_parse(file_name, map_style):
    def check_integer(integer, default):
        if integer not in option_dict:
            warnings.warn('%s is not set, defaults to %i' % (integer, default))
            option_dict[integer] = default
        elif not option_dict[integer].isdigit():
            line_split = option_dict[integer].split(',')
            if all(value.isdigit() for value in line_split):
                option_dict[integer] = [int(value) for value in line_split]
            else:
                warnings.warn('The %s given is not an integer, defaults to %i'
                              % (integer, default))
                option_dict[integer] = default
        else:
            option_dict[integer] = int(option_dict[integer])

    def check_pat(pat):
        if pat not in option_dict:
            raise NameError("Could not find any specification of %s\n" % pat)
        elif (re.search('[^ACGT]', option_dict[pat]) is not None and
              option_dict[pat] != 'NA'):
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

    for integer, default in [('norm_index_length', 0),
                             ('exp_index_length', 0),
                             ('map_fwd_index_length', 0),
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


def top_map(map_dict_in, max_dist):

    def sort_keys(key, map_dict_in):
        return (map_dict_in[key][0])

    map_dict_out = {}
    for bc in map_dict_in:
        this_map_dict = map_dict_in[bc]
        if len(this_map_dict) == 1:
            key = this_map_dict.keys()[0]
            reference_name, ori, start_pos = key
            total_reads, mapq_sum = this_map_dict[key]
            av_mapq = mapq_sum / total_reads
            freq1 = 1.0
            freq2 = 0.0
        else:
            total_reads = sum(value[0] for value in this_map_dict.values())
            sorted_key_list = sorted(this_map_dict.keys(),
                                     key=lambda elem: sort_keys(elem,
                                                                this_map_dict),
                                     reverse=True)
            top_key, \
                top_reads, \
                top_mapq, \
                sorted_key_list = refine_map(this_map_dict, sorted_key_list,
                                             max_dist)
            reference_name, ori, start_pos = top_key
            av_mapq = top_mapq / top_reads
            freq1 = top_reads / total_reads
            if len(sorted_key_list) > 0:
                second_key, \
                    second_reads, \
                    second_mapq, \
                    sorted_key_list = refine_map(this_map_dict,
                                                 sorted_key_list,
                                                 max_dist)
                freq2 = second_reads / total_reads
            else:
                freq2 = 0
        map_dict_out[bc] = [reference_name, ori, start_pos, total_reads,
                            av_mapq, freq1, freq2]
    return(map_dict_out)


def refine_map(this_map_dict, sorted_key_list, max_dist):
    top_key = sorted_key_list.pop(0)
    this_reads = this_map_dict[top_key][0]
    this_mapq = this_map_dict[top_key][1]
    i = 0
    while i < len(sorted_key_list):
        if (sorted_key_list[i][0] == top_key[0] and
                sorted_key_list[i][1] == top_key[1] and
                abs(sorted_key_list[i][2] - top_key[2]) < max_dist):
            close_key = sorted_key_list.pop(i)
            this_reads += this_map_dict[close_key][0]
            this_mapq += this_map_dict[close_key][1]
        i += 1
    return(top_key, this_reads, this_mapq, sorted_key_list)


def parse_sam(sam_file, max_soft_clip=5, min_first_match=10,
              remap_soft_clip=17):
    remap_list = []
    print sam_file
    map_dict = {}
    no_map_dict = {}
    for line in pysam.AlignmentFile(sam_file):
        name_split = line.query_name.split('_')
        if len(name_split) > 1:
            bc_this = name_split[1]
            if line.is_reverse:
                start_pos = line.reference_end
                ori = '-'
            else:
                start_pos = line.reference_start
                ori = '+'
            add_mapping = False
            if line.cigarstring is None:
                mapping = ('*', ori, 0)
                add_mapping = True
            else:
                if line.is_reverse:
                    first_cigar = line.cigartuples[-1]
                else:
                    first_cigar = line.cigartuples[0]
                mapping = (line.reference_name, ori, start_pos)
                if first_cigar[0] == 4:
                    if first_cigar[1] <= max_soft_clip:
                        add_mapping = True
                    if first_cigar[1] >= remap_soft_clip:
                        quality_dict = {'phred_quality':
                                        line.query_qualities[:first_cigar[1]]}
                        query_seq = line.query_sequence[:first_cigar[1]]
                        seq = SeqIO.SeqRecord(query_seq, line.query_name,
                                              description=line.query_name,
                                              letter_annotations=quality_dict)
                        remap_list.append(seq)
                if first_cigar[0] == 0 and first_cigar[1] >= 10:
                    add_mapping = True
            mapping_quality = line.mapping_quality
            if add_mapping:
                if bc_this in map_dict:
                    if mapping in map_dict[bc_this]:
                        map_dict[bc_this][mapping][0] += 1
                        map_dict[bc_this][mapping][1] += mapping_quality
                    else:
                        map_dict[bc_this][mapping] = [1, mapping_quality]
                else:
                    map_dict[bc_this] = {mapping: [1, mapping_quality]}
            else:
                if bc_this in no_map_dict:
                    if mapping in no_map_dict[bc_this]:
                        no_map_dict[bc_this][mapping][0] += 1
                        no_map_dict[bc_this][mapping][1] += mapping_quality
                    else:
                        no_map_dict[bc_this][mapping] = [1, mapping_quality]
                else:
                    no_map_dict[bc_this] = {mapping: [1, mapping_quality]}
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
    cmd_parser.add_argument('-u', '--useMagic', action='store_true',
                            dest='useMagic',
                            help=("use [] fields to specify multiple files. "
                                  "e.g.: 'TRIP_normalization_[12-22,24].fq' "
                                  "to specify all files numbered 12 to 22 "
                                  "and 24"))
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
    cmd_parser.add_argument('-p', '--count-partial', action='store_true',
                            dest='count_mutated',
                            help=('add the count of "mutated" barcodes to '
                                  'the count of the "genuine" barcode it was '
                                  'mutated from. "Genuine" barcodes are based '
                                  'on either bcFile or normalization reads'))
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


def run_top_map(map_dict_in_list, max_dist_list, cores):
    if cores > 1:
        pool = Pool(processes=cores)
        it = itertools.izip(itertools.repeat(top_map),
                            map_dict_in_list, max_dist_list)
        out = pool.map(function_star, it)
    else:
        out = []
        for i in range(0, len(map_dict_in_list)):
            out.append(top_map(map_dict_in_list[i], max_dist_list[i]))
    return out


def parse_bc_file(barcode_file):
    regex = re.compile(r' |\t')
    with open(barcode_file) as bc_file:
        for line in bc_file.readlines():
            yield re.split(regex, line.strip())[0:2]


def run_starcode(bc_dict, lev_dist, file_subset=[], bc_set=None, cores=1,
                 starcode_out=None):
    star_in_list = []
    if bc_set is not None:
        for bc in bc_set:
            star_in_list.append('%s\t%i' % (bc, 1000))
        for bc in bc_dict:
            if bc not in bc_set:
                star_in_list.append('%s\t%i' % (bc, 1))
    else:
        bc_set = set()
        for bc in bc_dict:
            if len(file_subset) > 0:
                value_list = [bc_dict[bc][f] for f in
                              file_subset if f in bc_dict[bc]]
            else:
                value_list = [value for value in bc_dict[bc].values()]
            norm_sum = sum(value_list)
            bc_set.add(bc)
            star_in_list.append('%s\t%i' % (bc, norm_sum))
    star_in = '\n'.join(star_in_list)
    with open('starcode_in.txt', 'w') as f:
        f.write(star_in)
    cmd_list = ['/home/NFS/users/c.leemans/Programs/starcode/starcode',
                '-t%i' % (cores), '-d%i' % lev_dist, '-s', '--print-clusters']
    pipe = subprocess.Popen(cmd_list, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE)
    stdout, stderr = pipe.communicate(input=star_in)
    if starcode_out is not None:
        with open(starcode_out, 'w') as f:
            f.write(stdout)
    for line in stdout.split('\n')[:-1]:
        line_list = line.split('\t')
        genuine_bc = line_list[0]
        if len(line_list) == 3:
            mutant_list = line_list[2].split(',')
            mutant_list = [bc for bc in mutant_list if bc != genuine_bc]
        else:
            mutant_list = []
        yield(genuine_bc, mutant_list)


def filter_barcodes(bc_dict, lev_dist, min_counts, norm_file_list=[],
                    bc_set=None, cores=1):
    exact_match_dict = {}
    part_match_dict = {}
    print datetime.now()
    star_in_list = []
    if bc_set is not None:
        for bc in bc_set:
            star_in_list.append('%s\t%i' % (bc, 1000))
        for bc in bc_dict:
            if bc not in bc_set:
                star_in_list.append('%s\t%i' % (bc, 1))
    else:
        bc_set = set()
        i = 0
        for bc in bc_dict:
            if len(norm_file_list) > 0:
                norm_base_list = (os.path.basename(norm_file)
                                  for norm_file in norm_file_list)
                norm_list = [bc_dict[bc][norm_base] for norm_base in
                             norm_base_list if norm_base in bc_dict[bc]]
            else:
                norm_list = [value for value in bc_dict[bc].values()]
            if sum(norm_list) > min_counts:
                norm_sum = sum(norm_list)
                bc_set.add(bc)
            else:
                norm_sum = 1
            if i < 20:
                print('%s\t%i' % (bc, norm_sum))
                i += 1
            star_in_list.append('%s\t%i' % (bc, norm_sum))

    star_in = '\n'.join(star_in_list)
    with open('starcode_in.txt', 'w') as f:
        f.write(star_in)
    cmd_list = ['/home/NFS/users/c.leemans/Programs/starcode/starcode',
                '-t%i' % (cores), '-d%i' % lev_dist, '-s', '--print-clusters']
    pipe = subprocess.Popen(cmd_list, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE)
    stdout, stderr = pipe.communicate(input=star_in)
    new_bc_dict = {}
    fake_set = set()
    for line in stdout.split('\n')[:-1]:
        line_list = line.split('\t')
        genuine_bc = line_list[0]
        if genuine_bc in bc_set:
            if genuine_bc in bc_dict:
                for f_name in bc_dict[genuine_bc]:
                    if f_name in exact_match_dict:
                        exact_match_dict[f_name] += bc_dict[genuine_bc][f_name]
                    else:
                        exact_match_dict[f_name] = bc_dict[genuine_bc][f_name]
                this_dict = {f_name: [count, 0] for f_name, count in
                             bc_dict[genuine_bc].items()}
            else:
                this_dict = {}
            if len(line_list) == 3:
                for fake_bc in line_list[2].split(','):
                    if fake_bc != genuine_bc:
                        if fake_bc in fake_set:
                            print fake_bc
                        else:
                            fake_set.add(fake_bc)
                        fake_bc_dict = bc_dict[fake_bc]
                        for f_name in fake_bc_dict:
                            if f_name in this_dict:
                                this_dict[f_name][1] += fake_bc_dict[f_name]
                            else:
                                this_dict[f_name] = [0, fake_bc_dict[f_name]]
                            if f_name in part_match_dict:
                                part_match_dict[f_name] += fake_bc_dict[f_name]
                            else:
                                part_match_dict[f_name] = fake_bc_dict[f_name]
            new_bc_dict[genuine_bc] = this_dict
    return new_bc_dict, exact_match_dict, part_match_dict, stdout


def function_star(args):
    return(args[0](*args[1:]))


class parsing_runner:
    def __init__(self, bc_len, out_dir, cores):
        self.bc_len = bc_len
        self.out_dir = out_dir
        self.it = []
        self.out = None
        self.pool = Pool(processes=cores)
        self.exp_it_index = None
        self.map_it_index = None
        self.map_structure_dict = {}
        self.norm_exp_structure_dict = {}

    def add_mapping(self, map_pair_list, index_len_list, map_pat1, map_pat2,
                    map_pat_rev):
        start_map = len(self.it)
        for index_len in index_len_list:
            if index_len not in self.map_structure_dict:
                structure_file = '%s/map_structure_index%i.txt' \
                                 % (self.out_dir, index_len)
                build_mapping_structure(index_len, self.bc_len, map_pat1,
                                        map_pat2, map_pat_rev, structure_file)
                self.map_structure_dict[index_len] = structure_file
        self.it.extend([(parse_reads, map_pair,
                         self.map_structure_dict[index_len], self.out_dir)
                       for map_pair, index_len
                       in zip(map_pair_list, index_len_list)])
        self.map_it_index = (start_map, len(self.it))

    def add_norm_exp(self, fastq_list, index_len_list, pat1, pat2):
        start_exp = len(self.it)
        for index_len in index_len_list:
            if index_len not in self.norm_exp_structure_dict:
                structure_file = '%s/norm_exp_structure_index%i.txt' \
                                 % (self.out_dir, index_len)
                build_exp_structure(index_len, self.bc_len, pat1, pat2,
                                    structure_file)
                self.norm_exp_structure_dict[index_len] = structure_file
        self.it.extend([(parse_reads, fastq,
                         self.norm_exp_structure_dict[index_len], self.out_dir)
                       for fastq, index_len
                       in zip(fastq_list, index_len_list)])
        self.exp_it_index = (start_exp, len(self.it))

    def run(self):
        out = self.pool.map(function_star, self.it)
        if self.map_it_index is not None:
            map_out = itertools.islice(out, self.map_it_index[0],
                                       self.map_it_index[1])
        else:
            map_out = None
        if self.exp_it_index is not None:
            exp_out = itertools.islice(out, self.exp_it_index[0],
                                       self.exp_it_index[1])
        else:
            exp_out = None
        return(map_out, exp_out)


def run_parse_sam(sam_list, cores, max_soft_clip=5, min_first_match=10,
                  remap_soft_clip=17):
    if cores > 1:
        pool = Pool(processes=cores)
        it = itertools.izip(itertools.repeat(parse_sam),
                            sam_list,
                            itertools.repeat(max_soft_clip),
                            itertools.repeat(min_first_match),
                            itertools.repeat(remap_soft_clip))
        out = pool.map(function_star, it)
    else:
        out = []
        for sam_file in sam_list:
            out.append(parse_sam(sam_file, max_soft_clip, min_first_match,
                                 remap_soft_clip))
    return out


def use_magic(file_string):
    sep_regex = r',(?![0-9,\-:]+\])'
    bracket_regex = re.compile(r'\[([0-9,:-]+)\]')
    split_regex = re.compile(r'[,\-:]')
    start_regex = re.compile(r'[0-9]+\Z')
    end_regex = re.compile(r'^[0-9]+')
    file_pattern_list = []
    start = 0
    for match in re.finditer(sep_regex, file_string):
        file_pattern_list.append(file_string[start:match.start()])
        start = match.end()
    file_pattern_list.append(file_string[start:])
    file_list = []
    for file_pattern in file_pattern_list:
        base_name = os.path.basename(file_pattern)
        bracket_match = re.search(bracket_regex, base_name)
        if bracket_match is not None:
            file_string = re.sub(bracket_regex, '%s', file_pattern)
            range_string = bracket_match.group(1)
            first_int = re.findall(end_regex, range_string)[0]
            range_list = [int(first_int)]
            for split_match in re.finditer(split_regex, range_string):
                split = split_match.group(0)
                start = re.findall(start_regex,
                                   range_string[:split_match.start()])[0]
                end = re.findall(end_regex,
                                 range_string[split_match.end():])[0]
                if split == ',':
                    range_list.append(int(end))
                elif split == '-' or split == ':':
                    range_list.extend(range(int(start) + 1, int(end) + 1))
            for num in range_list:
                file_name = file_string % num
                if os.path.exists(file_name):
                    file_list.append(file_name)
        else:
            file_list.append(file_pattern)
    return file_list


def combine_parse_out(parse_out, read_type):
    assert read_type in ('map', 'norm_exp')
    total_reads_dict = {}
    total_match_dict = {}
    bc_dict = {}
    # to keep order
    base_list = []
    for bc_count_dict, stats, file_base in parse_out:
        base_list.append(file_base)
        total_reads_dict[file_base] = stats['reads']
        if read_type == 'map':
            total_match_dict[file_base] = stats['reads_written']
        elif read_type == 'norm_exp':
            total_match_dict[file_base] = stats['pat2']
        for bc in bc_count_dict:
            if bc in bc_dict:
                bc_dict[bc][file_base] = bc_count_dict[bc]
            else:
                bc_dict[bc] = {file_base: bc_count_dict[bc]}
    return (bc_dict, (total_reads_dict, total_match_dict), base_list)


def starcode(bc_dict, starcode_out_name, lev_dist, min_counts,
             cores, read_type, bc_set=None):
    if bc_set is not None:
        filter_out = filter_barcodes(bc_dict, lev_dist, min_counts,
                                     bc_set=bc_set, cores=cores)
    else:
        filter_out = filter_barcodes(bc_dict, lev_dist, min_counts,
                                     norm_file_list, cores=cores)
    bc_dict, exact_match_dict, part_match_dict, stdout = filter_out
    with open('/'.join((out_dir, starcode_out_name)), 'w') as f_out:
        f_out.write(stdout)
    return(bc_dict, exact_match_dict, part_match_dict)


def get_base_noext(file_name):
    base = os.path.basename(file_name)
    return(os.path.splitext(base)[0])


def norm_and_exp(parse_out, norm_file_list, bc_file, out_dir, min_counts,
                 lev_dist, count_mutated, cores, verbose):
    if bc_file is not None:
        bc_set = set(bc for index, bc in parse_bc_file(bc_file))
    else:
        bc_set = None
    bc_dict, stats, base_list = combine_parse_out(parse_out, 'norm_exp')
    norm_base_list = [get_base_noext(file) for file in norm_file_list]
    starcode_out = '/'.join((out_dir, 'norm_exp_starcode.txt'))
    exact_match_dict = {base: 0 for base in base_list}
    part_match_dict = {base: 0 for base in base_list}
    filtered_bc_dict = {}
    for genuine_bc, mutant_list in run_starcode(bc_dict, lev_dist,
                                                norm_base_list, bc_set,
                                                cores, starcode_out):
        count_dict = bc_dict[genuine_bc]
        for base in count_dict:
            exact_match_dict[base] += count_dict[base]
        for mutant_bc in mutant_list:
            mutant_count = bc_dict[mutant_bc]
            for base_name in mutant_count:
                part_match_dict[base_name] += mutant_count[base_name]
                if count_mutated:
                    if base_name in count_dict:
                        count_dict[base_name] += mutant_count[base_name]
                    else:
                        count_dict[base_name] = mutant_count[base_name]
        filtered_bc_dict[genuine_bc] = count_dict
    with open('/'.join((out_dir, 'bc_count.txt')), 'w') as f_out:
        first_line_list = ['bc']
        first_line_list.extend(base_list)
        f_out.write('\t'.join(first_line_list))
        f_out.write('\n')
        for bc in filtered_bc_dict:
            line_list = [bc]
            this_count_dict = filtered_bc_dict[bc]
            for file_name in base_list:
                if file_name in this_count_dict:
                    count = this_count_dict[file_name]
                    line_list.append(str(count))
                else:
                    line_list.append('0')
            f_out.write('\t'.join(line_list))
            f_out.write('\n')
    write_stats('/'.join((out_dir, 'norm_exp_stats.txt')),
                stats[0], stats[1], exact_match_dict,
                part_match_dict, base_list)
    return bc_dict


def write_stats(stats_file, total_reads_dict, total_match_dict,
                exact_match_dict, part_match_dict, base_list):
    with open(stats_file, 'w') as f_out:
        first_line_list = ['stat']
        first_line_list.extend(base_list)
        f_out.write('\t'.join(first_line_list))
        f_out.write('\n')
        for stat, match_dict in (('total reads', total_reads_dict),
                                 ('barcodes', total_match_dict),
                                 ('exact', exact_match_dict),
                                 ('partly', part_match_dict)):
            line_list = [stat]
            print(match_dict)
            print(stat)
            for file_name in base_list:
                line_list.append(str(match_dict[file_name]))
            f_out.write('\t'.join(line_list))
            f_out.write('\n')


def filter_map_reads(file_base, out_dir, mutant_dict):
    fwd_in = gzip.open('%s/%s_1.fastq.gz' % (out_dir, file_base))
    rev_in = gzip.open('%s/%s_2.fastq.gz' % (out_dir, file_base))
    fwd_out_name = '%s/%s_fwd.fastq.gz' % (out_dir, file_base)
    rev_out_name = '%s/%s_rev.fastq.gz' % (out_dir, file_base)
    fwd_out = gzip.open(fwd_out_name, 'w')
    rev_out = gzip.open(rev_out_name, 'w')
    id_dict = {}
    for read_id, bc in parse_barcodes('%s/%s.barcode.txt.gz' %
                                      (out_dir, file_base)):
        read_id = read_id.split(' ')[0]
        if bc in mutant_dict:
            if count_mutated:
                id_dict[read_id] = mutant_dict[bc]
        else:
            id_dict[read_id] = bc
    fwd_seqIO = SeqIO.parse(fwd_in, 'fastq')
    rev_seqIO = SeqIO.parse(rev_in, 'fastq')
    for fwd_record, rev_record in itertools.izip(fwd_seqIO, rev_seqIO):
        if (rev_record.id in id_dict and
                len(fwd_record.seq) > 2 and
                len(rev_record.seq) > 2):
            barcode = id_dict[rev_record.id]
            old_id = fwd_record.id
            new_id = '_'.join((old_id, barcode))
            fwd_record.id = rev_record.id = new_id
            fwd_record.description = \
                fwd_record.description.replace(old_id, new_id)
            rev_record.description = \
                rev_record.description.replace(old_id, new_id)
            SeqIO.write(fwd_record, fwd_out, 'fastq')
            SeqIO.write(rev_record, rev_out, 'fastq')
    fwd_in.close()
    rev_in.close()
    fwd_out.close()
    rev_out.close()
    return(fwd_out_name, rev_out_name)


def run_filter_reads(base_list, out_dir, mutant_dict, cores):
    if cores > 1:
        pool = Pool(processes=cores)
        it = itertools.izip(itertools.repeat(filter_map_reads),
                            base_list,
                            itertools.repeat(out_dir),
                            itertools.repeat(mutant_dict))
        out = pool.map(function_star, it)
    else:
        out = []
        for file_base in base_list:
            out.append(filter_map_reads(file_base, out_dir, mutant_dict))
    return out


def map_fwd_rev(parse_out, exp_bc_dict, max_dist_for, max_dist_rev,
                lev_dist, count_mutated, bowtie_base, out_dir, cores):
    map_bc_dict, stats, base_list = combine_parse_out(parse_out, 'map')
    starcode_out = '/'.join((out_dir, 'mapping_starcode.txt'))
    exact_match_dict = {base: 0 for base in base_list}
    part_match_dict = {base: 0 for base in base_list}
    mutant_dict = {}
    print(datetime.now())
    for genuine_bc, mutant_list in run_starcode(map_bc_dict, lev_dist,
                                                bc_set=exp_bc_dict,
                                                cores=cores,
                                                starcode_out=starcode_out):
        if genuine_bc in map_bc_dict:
            count_dict = map_bc_dict[genuine_bc]
            for base_name in count_dict:
                exact_match_dict[base_name] += count_dict[base_name]
        for mutant_bc in mutant_list:
            if mutant_bc in map_bc_dict:
                mutant_count = map_bc_dict[mutant_bc]
                for base_name in mutant_count:
                    part_match_dict[base_name] += mutant_count[base_name]
                mutant_dict[mutant_bc] = genuine_bc

    fastq_list = run_filter_reads(base_list, out_dir, mutant_dict, cores)
    fwd_out_name = ','.join((fwd for fwd, rev in fastq_list))
    rev_out_name = ','.join((rev for fwd, rev in fastq_list))
    print(datetime.now())
    map_sam_list = ['/'.join([out_dir, 'samFor.sam']),
                    '/'.join([out_dir, 'samRev.sam'])]
    align_command = ("bowtie2 -p %i -t --very-sensitive -x %s -U %s -S %s"
                     % (cores, bowtie_base, fwd_out_name, map_sam_list[0]))
    os.system(align_command)

    align_command = ("bowtie2 -p %i -t --very-sensitive-local -x %s -U %s "
                     "-S %s" % (cores, bowtie_base, rev_out_name,
                                map_sam_list[1]))
    os.system(align_command)

    # map_fwd_out = run_mapping_fwd(map_fwd_list, out_dir, map_pat1, map_pat2,
    #                               bc_len, index_len, bc_dict, cores)
    # id_dict_list = [map_fwd_out[fwd_fastq][1][0] for fwd_fastq in map_fwd_list]

    # total_fwd_dict = {map_fwd: map_fwd_out[map_fwd][1][1]
    #                   for map_fwd in map_fwd_list}
    # hits_fwd_dict = {map_fwd: map_fwd_out[map_fwd][1][2]
    #                  for map_fwd in map_fwd_list}
    # genuine_fwd_dict = {map_fwd: map_fwd_out[map_fwd][1][3]
    #                     for map_fwd in map_fwd_list}
    # if verbose:
    #     print 'parsing reverse iPCR reads...'
    # map_rev_out = run_mapping_rev(map_rev_list, id_dict_list, map_pat_rev,
    #                               out_dir, cores)

    # hits_rev_dict = {map_rev: hits for map_rev, hits in
    #                  zip(map_rev_list, map_rev_out.values())}

    # map_fwd_string = ','.join(map_out for map_out, t in map_fwd_out.values())
    # map_rev_string = ','.join(map_rev_out.keys())

    print 'parsing sam'
    print map_sam_list
    sam_out = run_parse_sam(map_sam_list, cores)

    rev_map_dict, remap_list, no_map_dict = sam_out[1]
    if verbose:
        print 'found %i reads which need to be re-mapped' % len(remap_list)
    if len(remap_list) > 0:
        remap_fq = '/'.join([out_dir, 'remap.fq'])
        remap_sam = '/'.join([out_dir, 'samRev2.sam'])
        SeqIO.write(remap_list, remap_fq, 'fastq')
        align_command = ("bowtie2 -p %i -t --very-sensitive-local"
                         " -x %s -U %s -S %s" % (cores, bowtie_base, remap_fq,
                                                 remap_sam))
        os.system(align_command)
        remap_sam_out = parse_sam(remap_sam)
        for bc in remap_sam_out[0]:
            if bc in rev_map_dict:
                this_remap_dict = remap_sam_out[0][bc]
                this_map_dict = rev_map_dict[bc]
                for loc in this_remap_dict:
                    if loc in this_map_dict:
                        this_map_dict[loc] = [a + b for a, b in
                                              zip(this_remap_dict[loc],
                                                  this_map_dict[loc])]
                    else:
                        this_map_dict[loc] = this_remap_dict[loc]
            else:
                rev_map_dict[bc] = remap_sam_out[0][bc]
    fwd_map_dict = sam_out[0][0]
    top_map_out = run_top_map([fwd_map_dict, rev_map_dict],
                              [max_dist_for, max_dist_rev],
                              cores)

    header_tuple = ('chr', 'ori', 'pos', 't_reads', 'mapq', 'freq1', 'freq2')
    head_fwd = '\t'.join('%s_f' % head for head in header_tuple)
    head_rev = '\t'.join('%s_r' % head for head in header_tuple)
    unmapped_string = '\t'.join('NA' for i in header_tuple)
    if map_style == 'f':
        first_line = '\t'.join(('barcode', head_fwd))
    if map_style == 'r':
        first_line = '\t'.join(('barcode', head_rev))
    else:
        first_line = '\t'.join(('barcode', head_fwd, head_rev))
    line_template = '%s\t%s\t%i\t%i\t%g\t%g\t%g'
    with open('/'.join((out_dir, 'final_mapping.txt')), 'w') as file_out:
        file_out.write(first_line)
        file_out.write('\n')
        fwd_map_dict = top_map_out[0]
        rev_map_dict = top_map_out[1]
        for bc in fwd_map_dict:
            file_out.write(bc)
            if (map_style == 'f' or map_style == 'b'):
                file_out.write('\t')
                if (bc in fwd_map_dict):
                    fwd_map_list = fwd_map_dict[bc]
                    fwd_map_list[5] = round(fwd_map_list[5], 2)
                    fwd_map_list[6:] = [round(i, 2) for i in fwd_map_list[6:]]
                    file_out.write(line_template % tuple(fwd_map_list))
                else:
                    file_out.write(unmapped_string)
            if (map_style == 'r' or map_style == 'b'):
                file_out.write('\t')
                if (bc in rev_map_dict):
                    rev_map_list = rev_map_dict[bc]
                    rev_map_list[5] = round(rev_map_list[5], 2)
                    rev_map_list[6:] = [round(i, 3) for i in rev_map_list[6:]]
                    file_out.write(line_template % tuple(rev_map_list))
                else:
                    file_out.write(unmapped_string)
            file_out.write('\n')
    write_stats('/'.join((out_dir, 'mapping_stats.txt')),
                stats[0], stats[1], exact_match_dict,
                part_match_dict, base_list)


def extract_multi_hits(sam_in, fastq_out):
    seq_list = []
    for line in pysam.AlignmentFile(sam_in):
        if line.has_tag('XS'):
            quality_dict = {'phred_quality': line.query_qualities}
            seq = SeqIO.SeqRecord(line.query_sequence, line.query_name,
                                  description=line.query_name,
                                  letter_annotations=quality_dict)
            seq_list.append(seq)
    with open(fastq_out, 'w') as fout:
        SeqIO.write(seq_list, fout, 'fasta')


if __name__ == '__main__':
    options = get_arg_options()
    config_file = options.configuration_file
    if options.map_style is not None:
        map_style = options.map_style
    else:
        map_style = 'n'
    config_dict = config_parse(config_file, map_style)
    if options.normalization_file is not None:
        if options.useMagic:
            print options.normalization_file
            norm_file_list = use_magic(options.normalization_file)
        else:
            norm_file_list = options.normalization_file.split(',')
        norm_index_len = config_dict['norm_index_length']
        if type(norm_index_len) == int:
            norm_index_len_list = [norm_index_len for i in norm_file_list]
        else:
            norm_index_len_list = norm_index_len
    else:
        norm_file_list = []
        norm_index_len_list = []
    if options.expression_file is not None:
        if options.useMagic:
            print options.expression_file
            exp_file_list = use_magic(options.expression_file)
        else:
            exp_file_list = options.expression_file.split(',')
        exp_index_len = config_dict['exp_index_length']
        if type(exp_index_len) == int:
            exp_index_len_list = [exp_index_len for i in exp_file_list]
        else:
            exp_index_len_list = exp_index_len
    else:
        exp_file_list = []
        exp_index_len_list = []
    if options.mapping_forward is not None:
        if options.useMagic:
            map_fwd_list = use_magic(options.mapping_forward)
        else:
            map_fwd_list = options.mapping_forward.split(',')
        map_index_len = config_dict['map_fwd_index_length']
        if type(map_index_len) == int:
            map_index_len_list = [map_index_len for i in
                                  map_fwd_list]
        else:
            map_index_len_list = map_index_len
    else:
        map_fwd_list = []
    if options.mapping_reverse is not None:
        if options.useMagic:
            map_rev_list = use_magic(options.mapping_reverse)
        else:
            map_rev_list = options.mapping_reverse.split(',')
    else:
        map_rev_list = []
    count_mutated = options.count_mutated
    out_dir = options.output_directory
    bc_file = options.barcode_file
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

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
        debug_file.write('norm_index_length:\t%')
        debug_file.write(str(config_dict['norm_index_length']))
        debug_file.write('\n')
        debug_file.write('exp_index_length:\t%')
        debug_file.write(str(config_dict['exp_index_length']))
        debug_file.write('\n')
        debug_file.write('map_fwd_index_length:\t%')
        debug_file.write(str(config_dict['map_fwd_index_length']))
        debug_file.write('\n')
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
    runner = parsing_runner(config_dict['barcode_length'], out_dir,
                            config_dict['cores'])
    if (len(norm_file_list) + len(exp_file_list)) > 0:
        if verbose:
            print 'parsing normalization and expression reads'
        index_len_list = []
        index_len_list.extend(norm_index_len_list)
        index_len_list.extend(exp_index_len_list)
        pat1 = '' if config_dict['pat1'] == 'NA' else config_dict['pat1']
        # bc_dict = norm_and_exp(norm_file_list, exp_file_list, bc_file, out_dir,
        #                        config_dict['barcode_length'], pat1,
        #                        config_dict['pat2'], index_len_list,
        #                        config_dict['min_counts'],
        #                        config_dict['lev_dist'], count_mutated,
        #                        config_dict['cores'], verbose)

        runner.add_norm_exp(itertools.chain(norm_file_list, exp_file_list),
                            index_len_list, pat1, config_dict['pat2'])
        # build_exp_structure(config_dict['index_length'],
        #                     config_dict['barcode_length'], pat1,
        #                     config_dict['pat2'], 'tomtest/exp_structure.txt')
        # test = parse_reads(norm_file_list[0], 'tomtest/exp_structure.txt',
        #                    out_dir)
    else:
        if verbose:
            print 'parsing barcode file...'
        bc_dict = {key: {'index': value}
                   for value, key in parse_bc_file(bc_file)}
        # bc_tree = BKTree(Levenshtein.distance, bc_dict.keys())

    if len(map_fwd_list) > 0:
        if verbose:
            print 'parsing forward iPCR reads...'
        map_pat1 = '' if config_dict['map_pat1'] == 'NA' \
                   else config_dict['map_pat1']
        runner.add_mapping(zip(map_fwd_list, map_rev_list), map_index_len_list,
                           map_pat1, config_dict['map_pat2'],
                           config_dict['map_pat_rev'])
        # build_mapping_structure(config_dict['index_length'],
        #                         config_dict['barcode_length'],
        #                         map_pat1,
        #                         config_dict['map_pat2'],
        #                         config_dict['map_pat_rev'],
        #                         'tomtest/map_struc.txt')
        # test = parse_reads((map_fwd_list[0], map_rev_list[0]),
        #                    'tomtest/map_struc.txt', out_dir)
        # map_fwd_rev(map_fwd_list, map_rev_list, bc_dict,
        #             map_pat1, config_dict['map_pat2'],
        #             config_dict['map_pat_rev'], config_dict['barcode_length'],
        #             map_fwd_index_len_list, config_dict['max_dist_for'],
        #             config_dict['max_dist_rev'], config_dict['bowtie_base'],
        #             out_dir, config_dict['cores'])
    print 'parsing fastq files...'
    map_out, exp_out = runner.run()
    if exp_out is not None:
        print 'analysing normalization and expression reads...'
        bc_dict = norm_and_exp(exp_out, norm_file_list, bc_file, out_dir,
                               config_dict['min_counts'],
                               config_dict['lev_dist'], count_mutated,
                               config_dict['cores'], verbose)
    if map_out is not None:
        print 'analysing mapping reads...'
        map_fwd_rev(map_out, bc_dict, config_dict['max_dist_for'],
                    config_dict['max_dist_rev'], config_dict['lev_dist'],
                    count_mutated, config_dict['bowtie_base'], out_dir,
                    config_dict['cores'])
    # ind_len = 10
    # bc_len = 16
    # min_count = 5
    # pat1 = 'GTCACAAGGGCCGGCCACAACTCGAG'
    # pat2 = 'TGATCCTGCAGTG'
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
    # bc_count_dict = bc_extract_exp1(fastq_name, bc_len, pat1, pat2,
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
