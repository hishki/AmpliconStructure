import pysam
from numpy import genfromtxt

import collections
import bisect

READ_PATH = 'data/merged.bam'
SEGMENT_PATH = 'data/segments.csv'
SEGMENT_IN_READ_POS_THRESHOLD = 10


def get_reads(read_path):
    samfile = pysam.AlignmentFile(read_path, 'rb')
    all_reads = []
    for read in samfile.fetch():
        ref_name = read.reference_name.split(':')[0]
        ref_start = int(read.reference_name.split(':')[1].split('-')[0])
        ref_pos_to_query_pos = {}
        for query_pos, ref_pos in read.get_aligned_pairs(True):
            query_pos = read.query_length - query_pos - 1 if read.is_reverse else query_pos
            ref_pos_to_query_pos[ref_start + ref_pos] = query_pos
        ref_pos_to_query_pos = collections.OrderedDict(ref_pos_to_query_pos)
        all_reads.append({
            's': read.reference_start + ref_start,
            'e': read.reference_end + ref_start,
            'sup': read.is_supplementary,
            'chr': ref_name,
            'name': read.qname,
            'rev': read.is_reverse,
            'sec': read.is_secondary,
            'ref_pos_to_query_pos': ref_pos_to_query_pos
        })
    return all_reads


def map_read_by_segment(seg_path, reads):
    segments = genfromtxt(seg_path, delimiter='\t', dtype=str)
    seg_to_read = {}
    for segment in segments:
        _, name, chr, start, end = segment
        start, end = int(start), int(end)
        seg_to_read[name] = []
        for read in reads:
            if ((read['s'] <= start + SEGMENT_IN_READ_POS_THRESHOLD <= read['e']) or
                (read['s'] <= end <= read['e'] + SEGMENT_IN_READ_POS_THRESHOLD)) and chr.lower()[0:5] == read['chr'].lower()[0:5]:
                if start < read['s']:
                    start = read['s']
                if end >= read['e']:
                    end = read['e']-1
                # if start in read['ref_pos_to_query_pos'].keys():
                #     segment_start_pos_in_read = read['ref_pos_to_query_pos'][start]
                # else:
                #     print('wtf!!!')
                keys = list(read['ref_pos_to_query_pos'].keys())
                start = keys[max(bisect.bisect_left(keys, start) - 1, 0)]
                segment_start_pos_in_read = read['ref_pos_to_query_pos'][start]
                # if end in read['ref_pos_to_query_pos'].keys():
                end = keys[bisect.bisect_left(keys, end) - 1]
                # print(start, end)
                segment_end_pos_in_read = read['ref_pos_to_query_pos'][end]
                # else:
                #     read['ref_pos_to_query_pos'].keys()
                #     segment_end_pos_in_read = read['ref_pos_to_query_pos'][end-1]
                #     print('wtf!!!')
                if read['rev']:
                    segment_start_pos_in_read, segment_end_pos_in_read = segment_end_pos_in_read, segment_start_pos_in_read

                seg_to_read[name].append({
                    'segment_start_pos_in_read': segment_start_pos_in_read,
                    'segment_end_pos_in_read': segment_end_pos_in_read,
                    'rev': read['rev'],
                    'name': read['name']
                })
    return seg_to_read


reads = get_reads(READ_PATH)
mapped_data = map_read_by_segment(SEGMENT_PATH, reads)


if __name__ == "__main__":
    for key in mapped_data:
        print (key, mapped_data[key])
# print(mapped_data)



