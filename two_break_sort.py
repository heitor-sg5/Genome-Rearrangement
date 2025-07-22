def chromosome_to_cycle(chromosome):
    nodes = [0] * (2 * len(chromosome))
    for j in range(len(chromosome)):
        i = chromosome[j]
        if i > 0:
            nodes[2 * j] = 2 * i - 1
            nodes[2 * j + 1] = 2 * i
        else:
            nodes[2 * j] = -2 * i
            nodes[2 * j + 1] = -2 * i - 1
    return nodes

def cycle_to_chromosome(nodes):
    chromosome = []
    for j in range(len(nodes) // 2):
        if nodes[2 * j] < nodes[2 * j + 1]:
            chromosome.append(nodes[2 * j + 1] // 2)
        else:
            chromosome.append(-nodes[2 * j] // 2)
    return chromosome

def colored_edges(genome):
    edges = set()
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        nodes.append(nodes[0])
        for j in range(len(chromosome)):
            edges.add((nodes[2 * j + 1], nodes[2 * j + 2]))
    return edges

def two_break_on_genome_graph(edges, i0, i1, j0, j1):
    edges.discard((i0, i1))
    edges.discard((i1, i0))
    edges.discard((j0, j1))
    edges.discard((j1, j0))
    edges.add((i0, j0))
    edges.add((i1, j1))
    return edges

def find_and_merge(elements):
    parent = {x: x for x in elements}
    rank = {x: 0 for x in elements}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def merge(x, y):
        x_root = find(x)
        y_root = find(y)
        if x_root == y_root:
            return
        if rank[x_root] > rank[y_root]:
            parent[y_root] = x_root
        else:
            parent[x_root] = y_root
            if rank[x_root] == rank[y_root]:
                rank[y_root] += 1

    return find, merge

def group_nodes(edges):
    elements = set()
    for a, b in edges:
        elements.update([a, b])
        elements.update([a + 1 if a % 2 else a - 1])
        elements.update([b + 1 if b % 2 else b - 1])

    find, merge = find_and_merge(elements)

    for a, b in edges:
        merge(a, b)
        merge(a, a + 1 if a % 2 else a - 1)
        merge(b, b + 1 if b % 2 else b - 1)

    nodes_id = {x: find(x) for x in elements}
    return nodes_id

def build_edge_dict(edges, nodes_id):
    edge_dict = dict()
    for e in edges:
        id = nodes_id[e[0]]
        if id not in edge_dict:
            edge_dict[id] = dict()
        edge_dict[id][e[0]] = e[1]
        edge_dict[id][e[1]] = e[0]
    return edge_dict

def two_break_on_genome(genome, i0, i1, j0, j1):
    edges = two_break_on_genome_graph(colored_edges(genome), i0, i1, j0, j1)
    nodes_id = group_nodes(edges)
    edge_dict = build_edge_dict(edges, nodes_id)
    nodes_dict = dict()
    for id, edge in edge_dict.items():
        nodes_dict[id] = []
        curr0 = list(edge)[0]
        while edge:
            nodes_dict[id].append(curr0)
            if curr0 % 2 == 1:
                curr1 = curr0 + 1
            else:
                curr1 = curr0 - 1
            nodes_dict[id].append(curr1)
            new_node = edge[curr1]
            del edge[curr0]
            del edge[curr1]
            curr0 = new_node
    new_genome = []
    for nodes in nodes_dict.values():
        new_genome.append(cycle_to_chromosome(nodes))
    new_genome.sort(key=lambda x: abs(x[0]))
    return new_genome

def edge_from_non_trivial_cycle(edges, red_edges, blue_edges, blocks):
    elements = set()
    for a, b in edges:
        elements.update([a, b])

    find, merge = find_and_merge(elements)

    for a, b in edges:
        merge(a, b)

    nodes_id = {}
    nodes_sets = set()
    for a, b in edges:
        root = find(a)
        nodes_id[a] = root
        nodes_id[b] = root
        nodes_sets.add(root)

    cycles = len(nodes_sets)
    has_non_trivial_cycle = cycles != blocks
    removed = []

    if has_non_trivial_cycle:
        edge = None
        edge_dict = {}
        red_edge_dict = {}

        for a, b in edges:
            cid = nodes_id[a]
            edge_dict.setdefault(cid, {})[a] = b
            edge_dict[cid][b] = a
            if (a, b) in red_edges or (b, a) in red_edges:
                red_edge_dict.setdefault(cid, {})[a] = b
                red_edge_dict[cid][b] = a
            if edge is None and len(edge_dict[cid]) > 2 and (a, b) in blue_edges:
                edge = (a, b)
                edge_id = cid

        removed.append((edge[0], red_edge_dict[edge_id][edge[0]]))
        removed.append((edge[1], red_edge_dict[edge_id][edge[1]]))

    return has_non_trivial_cycle, removed

def shortest_rearrangement_scenario(P, Q):
    blocks = sum(len(chrom) for chrom in P)
    result = [P]
    red_edges = colored_edges(P)
    blue_edges = colored_edges(Q)
    breakpoint_graph = red_edges.union(blue_edges)
    has_non_trivial_cycle, removed = edge_from_non_trivial_cycle(breakpoint_graph, red_edges, blue_edges, blocks)
    while has_non_trivial_cycle:
        red_edges = two_break_on_genome_graph(red_edges, removed[0][0], removed[0][1], removed[1][0], removed[1][1])
        breakpoint_graph = red_edges.union(blue_edges)
        P = two_break_on_genome(P, removed[0][0], removed[0][1], removed[1][0], removed[1][1])
        has_non_trivial_cycle, removed = edge_from_non_trivial_cycle(breakpoint_graph, red_edges, blue_edges, blocks)
        result.append(P)
    return result

def two_break_distance_and_sort(P, Q):
    steps = shortest_rearrangement_scenario(P, Q)
    distance = len(steps) - 1
    return distance, steps

P = [[1, 2, 3, 4], [5, 6, 7, 8]]
Q = [[1, -3, -2], [4, 5], [6, 7, 8]]

distance, steps = two_break_distance_and_sort(P, Q)

i = 0
for step in steps:
    i += 1
    print(f'Step {i}: ', ''.join(['[' + ' '.join(f"{'+' if x > 0 else ''}{x}" for x in chrom) + ']' for chrom in step]))

print("2-Break Distance:", distance)
