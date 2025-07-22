# Genome Reaarangement Algorithms

This repository contains implementations of core algorithms for pairwise genome rearrangement analysis. These tools allow detection and comparison of synteny blocks between genomes, estimation of reversal distances, and computation of 2-break rearrangement scenarios. These methods support structural variation analysis and evolutionary comparisons of genome architecture.

---

## 🧬 What is Genome Rearrangement?

Genome rearrangement refers to large-scale structural changes in genome organization, such as inversions, translocations, fusions, and fissions. These rearrangements play an important role in evolution and can disrupt gene order while preserving sequence content. 

To understand how two genomes differ in structure, we can represent them as permutations of conserved blocks (synteny blocks) and compute the minimal set of operations (like reversals or 2-breaks) to transform one into another. Algorithms in this repository help identify these blocks and compute rearrangement distances.

---

## 📁 Files in This Repository

- `construct_synteny_blocks.py`: Identifies shared k-mers between two genomes, clusters them into synteny blocks, computes signed permutations, and generates a dot plot.
- `breakpoint_reversal_sort.py`: Greedily sorts a signed permutation into the identity permutation using reversals that reduce breakpoints.
- `two_break_sort.py`: Computes the 2-break distance between two genomes and lists the rearrangement steps transforming one into the other.
---

## ⚙️ How to Use

### 1. Prepare Input

Each script operates on a different form of genome input:

- `construct_synteny_blocks.py` prompts the user to select two genome FASTA files. It extracts DNA sequences and finds shared k-mers (default k = 30).
- `breakpoint_reversal_sort.py` requires a signed permutation, P, in the form of a list (e.g. [1, -2, -3]). You can generate this with `construct_synteny_blocks.py`.
- `two_break_sort.py` operates on genomes represented as lists of chromosomes, each a list of signed integers (permutations). 

### 2. Run the Algorithms

Each script will produce:

- An identity permutation P and a signed permutation Q, along with a genome dot-plot (P, Q), if using `construct_synteny_blocks.py`.
- A computed distance between genomes P and Q, and the intermediate steps between them. This will be reversals for `breakpoint_reversal_sort.py` and two-break operations for `two_break_sort.py`.

---

#### Synteny Block Constructor 

  bash
```construct_synteny_blocks.py```

#### Unichromosomal Reversal Sorting with Breakpoints

  bash
```breakpoint_reversal_sort.py```

#### Multichromosomal Two-break Sorting

  bash
```two_break_sort.py```

Parameters such as k (k-mer length), max_distance, and min_size in `construct_synteny_blocks.py` can be modified directly in the code to suit genome scale and resolution.

---

## 🧠 Algorithm Overviews

### Synteny Block Constructor 

- Builds an index of all k-mers (and reverse complements) in the second genome and identifies matching k-mers, classifying them by orientation (+/-).
- Constructs a proximity graph and finds connected components as candidate synteny blocks.
- Computes signed permutations by ordering blocks by their genome coordinates and orientation.
- Plots a genome dot-plot of all shared k-mers by orientation.
- Time complexity: O(n + m + s²)

### Unichromosomal Reversal Sorting with Breakpoints

- Counts breakpoints where adjacent elements are not consecutive.
- Greedily selects the reversal that maximally reduces breakpoints at each step.
- Applies the reversal, prints the intermediate permutation, and repeats, stopping when no breakpoints remain, indicating the identity permutation has been reached.
- Time complexity: O(n^3)

### Multichromosomal Two-break Sorting

- Converts genomes into breakpoint graphs using red (P) and blue (Q) edges.
- Identifies non-trivial cycles (edges not forming isolated cycles) and applies 2-breaks that reduce the number of such cycles.
- Updates the genome after each 2-break and repeats until P = Q, printing the transformation at each step.
- Time complexity: O(n^2)

---

## 🧪 Example Output

- Identity and Signed Permutations:

  Genome 1: [1, 2, 3, 4, 5]
  
  Genome 2: [1, -3, -2, 4, 5]
  
- Reversal Distance and Breakpoints:

  Step 1: [+1 +7 -9 +11 +10 +3 -2 -6 +5 -4 -8] | Breakpoints: 11
  Step 2: [+1 +2 -3 -10 -11 +9 -7 -6 +5 -4 -8] | Breakpoints: 9
  Step 3: [+1 +2 -3 -10 -11 +9 -7 -6 -5 -4 -8] | Breakpoints: 7
  Step 4 [+1 +2 +3 -10 -11 +9 -7 -6 -5 -4 -8] | Breakpoints: 6
  Step 5: [+1 +2 +3 +4 +5 +6 +7 -9 +11 +10 -8] | Breakpoints: 5
  Step 6 [+1 +2 +3 +4 +5 +6 +7 -11 +9 +10 -8] | Breakpoints: 4
  Step 7 [+1 +2 +3 +4 +5 +6 +7 +8 -10 -9 +11] | Breakpoints: 2
  Step 8: [+1 +2 +3 +4 +5 +6 +7 +8 +9 +10 +11] | Breakpoints: 0
  Reversal distance: 7

---

## 👤 Author

Heitor Gelain do Nascimento
Email: heitorgelain@outlook.com
GitHub: @heitor-sg5

---

## 📚 References

Bioinformatics Algorithms: An Active Learning Approach (Chapter 6) by
Phillip Compeau & Pavel Pevzner
https://bioinformaticsalgorithms.com
