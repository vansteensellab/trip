#!/usr/bin/python
from collections import Counter
from bitarray import bitarray
import re
import operator
from Bio.Data import IUPACData


class BKTreeCountDict:
    def __init__(self, distfn, words, compressor=None, decompressfn=None,
                 **kwargs):
        """
        Create a BK-tree with word counts from the given distance function and
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
        self.compressor = compressor
        it = iter(words)
        if compressor is not None:
            root = compressor.to_binary(it.next())
        self.tree = (compressor(root), {})
        self.count = {self.tree, 1}
        for i in it:
            self._add_word(self.tree, i)

    def _add_word(self, parent, word):
        pword, children = parent
        if self.compressor is not None:
            pword = self.compressor.from_binary(pword)
        d = self.distfn(word, pword, **self.kwargs)
        if self.compressor is not None:
            word = self.compressor.to_binary(word)
        if d == 0:
            self.count[word] += 1
        elif d in children:
            self._add_word(children[d], word)
        else:
            children[d] = (word, {})
            self.count[word] = 1

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
            if self.compressor is not None:
                pword = self.compressor.from_binary(pword)
            d = self.distfn(word, pword, **self.kwargs)
            results = []
            if d <= n:
                results.append((d, pword, self[pword]))
            for i in range(d - n, d + n + 1):
                child = children.get(i)
                if child is not None:
                    if self.compressor is not None:
                        child = self.compressor.from_binary(child)
                    results.extend(rec(child))
            return results

        # sort by distance
        return sorted(rec(self.tree))

    def __getitem__(self, word):
        """
        Just like a dictionary return the value stored for a key.
        In this case, the key is the word and the value is the count.
        """
        if self.compressor is not None:
            word = self.compressor.to_binary(word)
        return (self.count[word])


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

    def _add_word(self, parent, word):
        pword, children = parent
        d = self.distfn(word, pword, **self.kwargs)
        if d in children:
            self._add_word(children[d], word)
        else:
            children[d] = (word, {})

    def pop(self, word):
        def rec_child(new_parent, old_parent):
            for cword, children in old_parent.items():
                new_parent._add_word(new_parent, cword)
                if children is not None:
                    rec_child(new_parent, children)

        def rec(parent):
            pword, children = parent
            d = self.distfn(word, pword, **self.kwargs)
            if d == 0:
                if children != {}:
                    new_parent = children.pop(0)
                    rec_child(new_parent, children)
                    parent[0] = new_parent[0]
                    parent[1] = new_parent[1]
                return pword
            child = children.get(d)
            if child is not None:
                return(rec(child))
        return rec(self.tree)

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


class DNACompressor:
    def __init__(self):
        self.base_string = IUPACData.unambiguous_dna_letters
        self.iupac_fwd_dict = IUPACData.ambiguous_dna_values
        self.iupac_rev_dict = {v: k for k, v in
                               self.iupac_fwd_dict.iteritems()}
        self.baseN = len(self.base_string)

    def from_binary(self, bitString):
        bitArray = bitarray()
        bitArray.frombytes(bitString)
        dnaString = ''
        for i in range(0, len(bitArray), self.baseN):
            this_bit = bitArray[i:(i + self.baseN)]
            this_char = ''.join(self.base_string[j]
                                for j in range(0, self.baseN) if this_bit[j])
            if len(this_char) > 1:
                this_char = self.iupac_rev_dict[this_char]
            dnaString = ''.join([dnaString, this_char])
        return(dnaString)

    def to_binary(self, dnaString):
        bitString = bitarray()
        for char in dnaString:
            baseBit = [0 for i in range(0, self.baseN)]
            if char in self.iupac_fwd_dict:
                charIList = [self.base_string.index(c)
                             for c in self.iupac_fwd_dict[char]]
                for i in charIList:
                    baseBit[i] = 1
            else:
                raise SyntaxError("%s is not a valid DNA base" % (char))
            bitString.extend(baseBit)
        return bitString.tobytes()


def bc_extract_exp1(fastq, bc_len, constant1='',
                    constant2='', ind_len=0):
    bc_iterator = bc_extract(fastq, bc_len, constant1, constant2, ind_len)
    return(Counter(bc_iterator))


def bc_extract(fastq_name, bc_len, constant1='',
               constant2='', ind_len=0):
    len1 = len(constant1)
    len2 = len(constant2)
    regex = make_regex('exp', bc_len, constant1, constant2, ind_len)
    with open(fastq_name) as fastq:
        for line in fastq.readlines()[1::4]:
            match = re.match(regex, line)
            bc_this = None
            if match is not None:
                bc_this = match.group(1)
            else:
                hdist_1 = hamming_dist(line[ind_len:len1], constant1)
                if hdist_1 <= (len1 / 7):
                    start2 = (ind_len + len1 + bc_len)
                    if hamming_dist(line[start2:len1],
                                    constant2) <= (len2 / 7):
                        bc_this = line[(ind_len + len1):start2]
                    elif hamming_dist(line[(start2 - 1):len1],
                                      constant2) <= (len2 / 7):
                        bc_this = line[(ind_len + len1):(start2 - 1)]
                    elif hamming_dist(line[(start2 + 1):len1],
                                      constant2) <= (len2 / 7):
                        bc_this = line[(ind_len + len1):(start2 + 1)]
            if bc_this is not None:
                yield(bc_this)


def hamming_dist(string1, string2):
    min_len = min(len(string1), len(string2))
    hdist = sum(ch1 != ch2 for ch1, ch2 in
                zip(string1[:min_len], string2[:min_len]))
    return(hdist)


def make_regex(kind, bc_len, constant1='', constant2='', ind_len=0):
    pattern = r'^.{%i,%i}%s(.{%i,%i})%s%s'
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

    if (kind == 'exp'):
        end = '.+'
    elif (kind == 'map'):
        end = '(.+)'
    return re.compile(pattern % (min1, max1, constant1, bc_len - 1,
                                 bc_len + 1, constant2, end))


if __name__ == '__main__':
    ind_len = 10
    bc_len = 16
    constant1 = 'GTCACAAGGGCCGGCCACAACTCGAG'
    constant2 = 'TGATCCTGCAGTG'
    trip_folder = '/home/cleemans/SURFdrive/TRIP'
    fastq_base = '3884_1_BarcodedPool_NoIndex_TRIP_K562_KRAB_13.fq'
    fastq_name = '/'.paste(trip_folder, 'raw_data', fastq_base)
    countDict = bc_extract_exp1(fastq_name, bc_len, constant1, constant2,
                                ind_len)
    tree = BKTree(hamming_dist, countDict.keys())


    def sort_keys(key, countDict):
        return countDict[key]
    sorted_count = sorted(countDict.keys(),
                          key=lambda elem: sort_keys(elem, countDict),
                          reverse=True)
    i = 0
    while i < len(sorted_count):
        if round(i / 10) == (i / 10):
            print i
            print len(sorted_count)
        dna = sorted_count[i]
        neighbours = tree.query(dna, 1)
        for neighbour in neighbours:
            other_dna = neighbour[1]
            if other_dna != dna and other_dna in sorted_count:
                index = sorted_count.index(other_dna)
                dropped = sorted_count.pop(index)
        i += 1
