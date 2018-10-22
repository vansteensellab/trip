#!/usr/bin/env python
import subprocess


regions_fwd = snakemake.input['regions_fwd'][0]
bam_fwd = snakemake.input['fwd']
regions_rev = snakemake.input['regions_rev'][0]
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

mapping_dict = {}
with open(regions_fwd) as f_regions:
    line = f_regions.readline()
    for line in f_regions.readlines():
        line_split = line.strip().split()
        mapping_dict[line_split[0]] = {'fwd':line_split[1:7]}

with open(regions_rev) as f_regions:
    line = f_regions.readline()
    for line in f_regions.readlines():
        line_split = line.strip().split()
        bc = line_split[0]
        if bc not in mapping_dict:
            mapping_dict[bc] = {'rev':line_split[1:7]}
        else:
            mapping_dict[bc]['rev'] = line_split[1:7]

key_list = list(genomes.keys())
s0 = key_list[0]
s1 = key_list[1]

with open(output, 'w') as f_out:

    print(("barcode\tchr_f\tori_f\tpos_f\tt_reads_f\tmapq_f\tfreq1_f\tfreq2_f"
           "\tchr_r\tori_r\tpos_r\tt_reads_r\tmapq_r\tfreq1_r\tfreq2_r"
           "\t%s_mut\t%s_mut\tallele") % (s0, s1), file=f_out)
    for bc in mapping_dict:
        command_list = []
        mutation_count = {}
        info_list = ['NA' for i in range(0,14)]
        info_list.insert(0, bc)
        for strand in mapping_dict[bc]:
            map_list = mapping_dict[bc][strand]
            if strand == 'fwd':
                info_list[1:8] = map_list
            else:
                info_list[8:15] = map_list
        for s in genomes:
            for strand in mapping_dict[bc]:
                map_list = mapping_dict[bc][strand]
                seqname = map_list[0]
                ori = map_list[1]
                pos = int(map_list[2])
                start = pos if ori=='+' else pos - 230
                end = pos + 230 if ori=='+' else pos
                cmd = pileup_format.format(samtools, seqname, start, end,
                                           genomes[s], bam_fwd)
                command_list.append("{} | awk '{}'".format(cmd, awk))
            out_list = [run_shell(cmd) for cmd in command_list]
            mutation_count[s] = sum(int(b) for b in out_list)

        if mutation_count[s0] < mutation_count[s1]:
            allele = s0
        elif mutation_count[s1] < mutation_count[s0]:
            allele = s1
        else:
            allele = 'ambiguous'

        info_list.extend((str(mutation_count[s0]),
                          str(mutation_count[s1]),
                          allele))
        print('\t'.join(info_list), file=f_out)
