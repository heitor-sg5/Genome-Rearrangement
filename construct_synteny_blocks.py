from collections import defaultdict
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from collections import defaultdict, deque

def pattern_to_number(pattern):
    num = 0
    for symbol in pattern:
        num = 4 * num + symbol_to_number(symbol)
    return num

def symbol_to_number(symbol):
    return {'A': 0, 'C': 1, 'G': 2, 'T': 3}[symbol]

def reverse_complement_number(kmer_num, k):
    rc = 0
    for _ in range(k):
        rc = (rc << 2) | (3 - (kmer_num & 3))
        kmer_num >>= 2
    return rc

def find_shared_kmers(text1, text2, k):
    index2 = defaultdict(list)
    for j in range(len(text2) - k + 1):
        kmer = text2[j:j+k]
        kmer_num = pattern_to_number(kmer)
        index2[kmer_num].append(j)

    shared = []
    for i in range(len(text1) - k + 1):
        kmer = text1[i:i+k]
        kmer_num = pattern_to_number(kmer)
        rc_kmer_num = reverse_complement_number(kmer_num, k)

        for j in index2.get(kmer_num, ()):
            shared.append((i, j, '+'))
        for j in index2.get(rc_kmer_num, ()):
            shared.append((i, j, '-'))

    return shared

def build_synteny_graph(shared_kmers, max_distance):
    bin_size = max_distance
    bins = defaultdict(list)
    node_positions = []

    for idx, (i, j, _) in enumerate(shared_kmers):
        bin_x = i // bin_size
        bin_y = j // bin_size
        bins[(bin_x, bin_y)].append(idx)
        node_positions.append((i, j))
    adj = defaultdict(list)

    for idx, (i, j, _) in enumerate(shared_kmers):
        bin_x = i // bin_size
        bin_y = j // bin_size
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                neighbor_bin = (bin_x + dx, bin_y + dy)
                for neighbor_idx in bins.get(neighbor_bin, []):
                    if neighbor_idx == idx:
                        continue
                    i2, j2 = node_positions[neighbor_idx]
                    if abs(i - i2) <= max_distance and abs(j - j2) <= max_distance:
                        adj[idx].append(neighbor_idx)
    return adj

def find_connected_components(adj):
    visited = set()
    components = []
    for node in adj:
        if node not in visited:
            comp = []
            queue = deque([node])
            visited.add(node)
            while queue:
                current = queue.popleft()
                comp.append(current)
                for neighbor in adj[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            components.append(comp)
    return components

def synteny_blocks(shared_kmers, max_distance, min_size):
    adj = build_synteny_graph(shared_kmers, max_distance)
    comps = find_connected_components(adj)
    blocks = []
    for comp in comps:
        if len(comp) >= min_size:
            block = [shared_kmers[idx] for idx in comp]
            blocks.append(block)
    return blocks

def signed_permutations(blocks):
    metas = []
    for block in blocks:
        is_ = [i for i, _, _ in block]
        js = [j for _, j, _ in block]
        orients = [o for *_, o in block]
        avg_i = sum(is_) / len(is_)
        avg_j = sum(js) / len(js)
        sign = '+' if orients.count('+') >= orients.count('-') else '-'
        metas.append({'avg_i': avg_i, 'avg_j': avg_j, 'sign': sign})

    metas_sorted_by_p = sorted(metas, key=lambda x: x['avg_i'])
    for idx, m in enumerate(metas_sorted_by_p, start=1):
        m['id'] = idx

    perm1 = [m['id'] for m in metas_sorted_by_p]

    metas_sorted_by_q = sorted(metas, key=lambda x: x['avg_j'])
    perm2 = [m['id'] if m['sign'] == '+' else -m['id'] for m in metas_sorted_by_q]

    return perm1, perm2

def plot_dotplot(shared_kmers, k, genome1_len, genome2_len):
    x_f, y_f = [], []
    x_r, y_r = [], []
    
    for i, j, orientation in shared_kmers:
        if orientation == '+':
            x_f.append(i / 1000)  
            y_f.append(j / 1000)
        else:
            x_r.append(i / 1000)
            y_r.append(j / 1000)
    
    plt.figure(figsize=(10, 8))
    plt.scatter(x_f, y_f, color='red', s=10, label='(+) Matches')
    plt.scatter(x_r, y_r, color='blue', s=10, label='(-) Matches')
    plt.xlabel('P (kb)')
    plt.ylabel('Q (kb)')
    plt.title(f'Shared {k}-mers Dot-Plot')
    plt.legend()
    plt.xlim([0, genome1_len / 1000])
    plt.ylim([0, genome2_len / 1000])
    plt.grid(True)
    plt.show(block=True)

def extract_fasta_sequence(filename):
    with open(filename, 'r') as file:
        return ''.join(line.strip() for line in file if not line.startswith(">"))

def select_and_process_files():
    root = tk.Tk()
    root.withdraw()
    
    print("Select the first genome FASTA file:")
    file_path1 = filedialog.askopenfilename(
        title="Select the FIRST FASTA file",
        filetypes=[
            ("FASTA/FNA files", "*.fasta *.fa *.fna *.txt"),
            ("All files", "*.*"),
        ],
    )
    
    if not file_path1:
        print("First file not selected.")
        return None, None

    print("Select the second genome FASTA file:")
    file_path2 = filedialog.askopenfilename(
        title="Select the second FASTA file",
        filetypes=[
            ("FASTA/FNA files", "*.fasta *.fa *.fna *.txt"),
            ("All files", "*.*"),
        ],
    )

    if not file_path2:
        print("Second file not selected.")
        return None, None

    text1 = extract_fasta_sequence(file_path1)
    text2 = extract_fasta_sequence(file_path2)

    print("Both genomes successfully loaded.")
    return text1, text2

text1, text2 = select_and_process_files()

if text1 and text2:
    k = 30
    shared_kmers = find_shared_kmers(text1, text2, k)

    blocks = synteny_blocks(shared_kmers, max_distance=2000, min_size=500)
    perm1, perm2 = signed_permutations(blocks)

    print(f"Signed Permutations:\nGenome 1: {perm1}\nGenome 2: {perm2}")
    
    plot_dotplot(shared_kmers, k, len(text1), len(text2))

#Genome 1: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87]
#Genome 2: [1, 2, 3, 4, 77, 71, 62, -45, -28, 58, 11, 12, 13, 14, 23, 16, 17, 18, -20, -19, 21, 22, 15, 24, 43, -8, -79, -60, -69, -53, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, -57, -81, 25, -66, -72, -6, 47, -74, 49, 50, 51, 52, 61, -30, 70, 78, 10, -44, 59, 68, 55, 7, 64, -26, 82, -42, 67, 54, 63, 9, -27, 80, -41, -48, 75, 76, 5, -29, 56, -46, 65, 73, 83, 84, 85, 86, 87]
