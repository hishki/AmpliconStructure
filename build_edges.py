from pacbio_read import mapped_data
from numpy import genfromtxt

BEGIN_BEFORE_END_THRESHOLD = 10
BEGIN_AFTER_END_THRESHOLD = 1500
DOT_FILENAME = 'data/dot.gv'
SHORT_READ_EDGES_FILENAME = 'data/short_read_edges.csv'
short_read_edges = genfromtxt(SHORT_READ_EDGES_FILENAME, delimiter=' ', dtype=str)


def is_edge_exist(v, u, dir_v, dir_u):
    for edges in short_read_edges:
        if min(edges[0], edges[1]) == min(dir_v+v, dir_u+u) and max(edges[0], edges[1]) == max(dir_v+v, dir_u+u):
            return True
    return False


def build_edges(begin_before_end_threshold, begin_after_end_threshold):
    edges = {}
    for seg_x_name, seg_x_val in mapped_data.items():
        edges[seg_x_name] = {}
        for seg_y_name, seg_y_val in mapped_data.items():
            for read_x in seg_x_val:
                for read_y in seg_y_val:
                    if read_x['name'] == read_y['name']:
                        if read_x['segment_end_pos_in_read'] - begin_before_end_threshold <\
                           read_y['segment_start_pos_in_read'] <\
                           read_x['segment_end_pos_in_read'] + begin_after_end_threshold:
                                x_dir = 'b' if read_x['rev'] else 'e'
                                y_dir = 'e' if read_y['rev'] else 'b'
                                # if is_edge_exist(seg_x_name, seg_y_name, x_dir, y_dir):
                                if seg_y_name < seg_x_name:
                                    edges[seg_y_name][seg_x_name] = (y_dir, x_dir)
                                else:
                                    edges[seg_x_name][seg_y_name] = (x_dir, y_dir)
    return edges


def compare_short_long_edges(short_edges, long_edges):
    set_short_edges = set()
    set_long_edges = set()
    for edge in short_edges:
        set_short_edges.add((min(edge[0], edge[1]), max(edge[0], edge[1])))

    for v, v_edges in long_edges.items():
        for u, dir in v_edges.items():
            x, y = dir[0] + v, dir[1] + u
            x, y = min(x, y), max(x, y)
            set_long_edges.add((x, y))

    print("***intersect***")
    print(set_short_edges.intersection(set_long_edges))
    print("***diff***")
    print(set_short_edges.difference(set_long_edges))


def convert_to_dot(edges, filename):
    with open(filename, 'w') as f:
        f.write('graph graphname {\n')
        for v, v_edges in edges.items():
            if len(v_edges) > 0:
                pass
                f.write('{}{} -- {}{} [color=blue];\n'.format('b', v, 'e', v))
            for u, dir in v_edges.items():
                if min(int(v), int(u)) <= 48 < max(int(v), int(u)):
                    f.write('{}{} -- {}{} [color=green];\n'.format(dir[0], v, dir[1], u))
                else:
                    f.write('{}{} -- {}{};\n'.format(dir[0], v, dir[1], u))
        f.write('}')



def get_opposite_dir(dir):
    if dir == 'b':
        return 'e'
    return 'b'


def dfs(v, dir, is_parent_v, edges, path):
    mark[v][dir][is_parent_v] = 1
    path.append((dir+v, is_parent_v))
    print(path)
    if not is_parent_v:
        if mark[v][get_opposite_dir(dir)][not is_parent_v] == 0:
            # if (get_opposite_dir(dir)+v, not is_parent_v) in path:
            dfs(v, get_opposite_dir(dir), not is_parent_v, edges, path)
        elif mark[v][get_opposite_dir(dir)][not is_parent_v] == 1:
            print(path, get_opposite_dir(dir)+v)
    else:
        for u in edges[v]:
            dir_v, dir_u = edges[v][u]
            if dir_v == dir and v != u:
                if mark[u][dir_u][not is_parent_v] == 0:
                    # if (dir_u+u, not is_parent_v) not in path:
                    dfs(u, dir_u, not is_parent_v, edges, path)
                elif mark[u][dir_u][not is_parent_v] == 1:
                    print(path, dir_u + u)
    path.pop()
    mark[v][dir][is_parent_v] = 2
    # print('after', v, dir, path)


edges = build_edges(BEGIN_BEFORE_END_THRESHOLD, BEGIN_AFTER_END_THRESHOLD)
mark = {v: {'b': [0, 0], 'e': [0, 0]} for v in edges.keys()}

for v in edges.keys():
    for _dir in ('b', 'e'):
        if not mark[v][_dir][0]:
            dfs(v, _dir, False, edges, [])
        print(v, _dir)

compare_short_long_edges(short_read_edges, edges)
convert_to_dot(edges, DOT_FILENAME)

# for key, value in edges.items():
#     print(key, value)
