from pacbio_read import mapped_data

BEGIN_BEFORE_END_THRESHOLD = 100
BEGIN_AFTER_END_THRESHOLD = 100

edges = {}

for seg_x_name, seg_x_val in mapped_data.items():
    edges[seg_x_name] = {}
    for seg_y_name, seg_y_val in mapped_data.items():
        for read_x in seg_x_val:
            for read_y in seg_y_val:
                if read_x['name'] == read_y['name']:
                    if read_x['segment_end_pos_in_read'] - BEGIN_BEFORE_END_THRESHOLD <\
                       read_y['segment_start_pos_in_read'] <\
                       read_x['segment_end_pos_in_read'] + BEGIN_AFTER_END_THRESHOLD:
                            x_dir = '-' if read_x['rev'] else '+'
                            y_dir = '-' if read_y['rev'] else '+'
                            edges[seg_x_name][seg_y_name] = (x_dir, y_dir)

for key, value in edges.items():
    print(key, value)
