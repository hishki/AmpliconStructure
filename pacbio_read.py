import pysam
from numpy import genfromtxt

samfile = pysam.AlignmentFile('data/merged.bam', 'rb')

first_read_pos = samfile.fetch().__next__().pos

read_pos = []


for read in samfile.fetch():
    ref_name = read.reference_name.split(':')[0]
    ref_start = int(read.reference_name.split(':')[1].split('-')[0])
    read_pos.append({'s': read.reference_start + ref_start, 'e': read.reference_end + ref_start, 'sup': read.is_supplementary,
                     'chr': ref_name, 'name': read.qname, 'seg': [], 'rev': read.is_reverse, 'sec': read.is_secondary,
                     })

print(read_pos)

segments = genfromtxt('data/segments.csv', delimiter='\t', dtype=str)

for segment in segments:
    _, name, chr, start, end = segment
    for read in read_pos:
        if read['s'] < int(start)+1000 and read['e'] > int(end)-1000 and chr.lower()[0:5] == read['chr'].lower()[0:5]:
            read['seg'].append(segment)


print(read_pos)