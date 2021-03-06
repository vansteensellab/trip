#!/usr/bin/env python
import subprocess


regions = snakemake.input['regions'][0]
bam_fwd = snakemake.input['fwd']
bam_rev = snakemake.input['rev']
genomes = snakemake.params['genomes']

output = snakemake.output[0]

samtools = snakemake.params['samtools']

# bam_fwd = 'cl20170619_tn5/results/sorted/CM1420_fwd.bam'
# bam_rev = 'cl20170619_tn5/results/sorted/CM1420_rev.bam'
# regions = 'cl20170619_tn5/results/regions/CM1420_combined.txt'
# genomes =  {'CAST': '/home/NFS/users/c.leemans/data/GRCm38/GRCm38_CAST.fa',
#             '129S1': '/home/NFS/users/c.leemans/data/GRCm38/GRCm38_129S1.fa'}
# samtools = '/home/NFS/users/l.pagie/usr/local/src/miniconda2/bin/samtools'


pileup_format = "{} mpileup -t AD -v -d 50 -u -r {}:{}-{} -f {} {}"
awk = ('BEGIN {m=0}'
       '{'
       '    if ($1 ~ /^[^#]/){'
       '        match($8, "QS=([0-9.]*)",arr); '
       '        if (arr[1] < 0.5){'
       '            m += 1'
       '        }'
       '    }'
       '}'
       'END {'
       '    printf "%i", m'
       '}')



def run_shell(command):
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        outs, errs = p.communicate(timeout=15)
    except subprocess.TimeoutExpired:
        p.kill()
        outs, errs = p.communicate()
    return(outs)

with open(regions) as f_regions:
    with open(output, 'w') as f_out:
        print('#seqnames\tstart\tend\tori\tallele\tCAST_mut\t129S1_mut\tdepth_fwd\tdepth_rev', file=f_out)
        for line in f_regions.readlines():
            line_split = line.strip().split()
            mutation_count = {}
            for s in genomes:
                cmd_fwd = pileup_format.format(samtools, line_split[0],
                                               line_split[1], line_split[2],
                                               genomes[s], bam_fwd)
                cmd_rev = pileup_format.format(samtools, line_split[5],
                                               line_split[6], line_split[7],
                                               genomes[s], bam_rev)
                command_list = ["{} | awk '{}'".format(sam, awk)
                                for sam in (cmd_fwd, cmd_rev)]
                out_list = [run_shell(cmd) for cmd in command_list]
                mutation_count[s] = sum(int(b) for b in out_list)
            if line_split[3] == '+':
                middle = (int(line_split[1]), int(line_split[7]))
            else:
                middle = (int(line_split[2]), int(line_split[6]))
            insertion_start = min(middle)
            insertion_end = max(middle)
            key_list = list(genomes.keys())
            s0 = key_list[0]
            s1 = key_list[1]
            if mutation_count[s0] < mutation_count[s1]:
                allele = s0
            elif mutation_count[s1] < mutation_count[s0]:
                allele = s1
            else:
                allele = 'ambiguous'
            print('\t'.join((line_split[0], str(insertion_start),
                             str(insertion_end), line_split[3], allele,
                             str(mutation_count[s0]),
                             str(mutation_count[s1]), line_split[4],
                             line_split[9])), file=f_out)
