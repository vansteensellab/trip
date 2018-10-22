import pysam
import random
import itertools
import math
import re


sam_list = [pysam.AlignmentFile(b, 'rb', check_sq=False) for b in snakemake.input.bam]
sam_comb = pysam.AlignmentFile(snakemake.output[0], "wb", template=sam_list[0])
for line_tuple in itertools.zip_longest(*sam_list):
    match_list = [-math.inf if line.is_unmapped else line.get_tag('AS')
                  for line in line_tuple]
    line_list = [line_tuple[i] for i in range(0,2) if match_list[i]==max(match_list)]
    if len(line_list) == 1:
        sam_comb.write(line_list[0])
    else:
        sam_comb.write(line_list[round(random.random())])
[b.close() for b in sam_list]
sam_comb.close()
