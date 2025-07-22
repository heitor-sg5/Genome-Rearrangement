def count_breakpoints(perm):
    return sum(1 for i in range(len(perm) - 1) if perm[i + 1] - perm[i] != 1)

def apply_reversal(perm, i, j):
    return perm[:i] + [-x for x in perm[i:j + 1][::-1]] + perm[j + 1:]

def format_perm(perm):
    return '[' + ' '.join(f"{'+' if x > 0 else ''}{x}" for x in perm) + ']'

def reversal_sort_with_breakpoints(p):
    n = len(p)
    q = [0] + p + [n + 1]
    reversals = 0

    print("Step 1: ", format_perm(q[1:-1]), "| Breakpoints:", count_breakpoints(q))

    while count_breakpoints(q) > 0:
        best_q = None
        best_breaks = count_breakpoints(q)

        for i in range(1, len(q) - 1):
            for j in range(i, len(q) - 1):
                new_q = apply_reversal(q, i, j)
                new_breaks = count_breakpoints(new_q)
                if new_breaks < best_breaks:
                    best_breaks = new_breaks
                    best_q = new_q
                    if best_breaks == count_breakpoints(q) - 2:
                        break
            if best_q and best_breaks == count_breakpoints(q) - 2:
                break

        if best_q:
            q = best_q
            reversals += 1
            print(f"Step {reversals + 1}: ", format_perm(q[1:-1]), "| Breakpoints:", count_breakpoints(q))
        else:
            raise RuntimeError("Stuck: no reversal reduces breakpoints.")

    return reversals

p = [1, 7, -9, 11, 10, 3, -2, -6, 5, -4, -8]
print(f"Reversal distance: {reversal_sort_with_breakpoints(p)}")
