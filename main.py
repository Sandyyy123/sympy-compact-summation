"""
Compact-summation prototype for the Upwork job.

Client's ask, decoded:
  - A summation Sum(f(k), (k, lo, hi)) whose term f(k) is non-zero only at a
    sparse, spread-out set of indices (the "stored values", e.g. produced via
    .subs on the term).
  - Goal: MOVE those stored values so they sit at CONSECUTIVE indices at the
    start of a target range (e.g. 1..1000): 1st stored value -> index 1,
    2nd -> index 2, ... until the stored values are exhausted; every remaining
    index in the range returns 0.
  - The result must still be a real SymPy object: Sum() over it evaluates, and
    term.subs(k, i) returns the packed value at index i.

This file proves the operation works end-to-end and stays symbolic.
"""
import sympy as sp
from sympy import Sum, Piecewise, KroneckerDelta, symbols, Idx, IndexedBase, S

k = symbols('k', integer=True)


# ----------------------------------------------------------------------
# 1. Extract the stored (non-zero) values from a sparse symbolic term,
#    in index order, over the original range.
# ----------------------------------------------------------------------
def stored_values(term, index, lo, hi):
    """Return [(orig_index, value), ...] for indices where term != 0."""
    out = []
    for i in range(int(lo), int(hi) + 1):
        v = sp.simplify(term.subs(index, i))
        if v != 0:
            out.append((i, v))
    return out


# ----------------------------------------------------------------------
# 2a. Build a COMPACTED term as a Piecewise: value j lives at index j.
# ----------------------------------------------------------------------
def compact_piecewise(values, index, start=1):
    """values = [v1, v2, ...]; place them at start, start+1, ... else 0."""
    pieces = [(v, sp.Eq(index, start + j)) for j, v in enumerate(values)]
    pieces.append((S.Zero, True))          # everything else -> 0
    return Piecewise(*pieces)


# ----------------------------------------------------------------------
# 2b. Same thing built with KroneckerDelta (pure algebra, no Piecewise) --
#     handy when the client wants to keep manipulating it symbolically.
# ----------------------------------------------------------------------
def compact_kronecker(values, index, start=1):
    return sum(v * KroneckerDelta(index, start + j) for j, v in enumerate(values))


# ----------------------------------------------------------------------
# 2c. Same thing as an IndexedBase table -- closest to the client's
#     "indexes 1..1000" mental model; unset indices read as 0.
# ----------------------------------------------------------------------
def compact_indexed(values, name='a', start=1, length=1000):
    a = IndexedBase(name)
    table = {start + j: v for j, v in enumerate(values)}
    def lookup(i):
        return table.get(int(i), S.Zero)
    return a, table, lookup


# ======================================================================
# DEMO
# ======================================================================
if __name__ == '__main__':
    # A term that is non-zero only at scattered indices 3, 17, 42, 99, 500
    # (built here with KroneckerDelta so it's a genuine sparse SymPy term;
    #  in the client's code these come from .subs on their summand).
    original = (7*KroneckerDelta(k, 3) + 4*KroneckerDelta(k, 17)
                + 9*KroneckerDelta(k, 42) - 2*KroneckerDelta(k, 99)
                + 5*KroneckerDelta(k, 500))

    lo, hi = 1, 1000
    stored = stored_values(original, k, lo, hi)
    print("Stored values found (scattered):")
    for i, v in stored:
        print(f"   original index {i:>4}  ->  value {v}")

    vals = [v for _, v in stored]

    # --- compact into consecutive indices 1,2,3,... ---
    pw = compact_piecewise(vals, k, start=1)
    kd = compact_kronecker(vals, k, start=1)
    a, table, lookup = compact_indexed(vals, start=1, length=1000)

    print("\nAfter compaction -- values now sit at consecutive indices:")
    for i in range(1, len(vals) + 3):          # show packed region + a couple zeros
        got = pw.subs(k, i)
        print(f"   index {i}: Piecewise.subs -> {str(got):>3} | "
              f"Kronecker.subs -> {str(kd.subs(k, i)):>3} | table -> {lookup(i)}")

    # --- the three representations must agree at every index, and Sum works ---
    N = 1000
    # KroneckerDelta form evaluates SYMBOLICALLY under Sum().doit();
    # Piecewise-with-Eq stays unevaluated in SymPy, so evaluate it explicitly.
    s_pw = sum(pw.subs(k, i) for i in range(1, N + 1))
    s_kd = Sum(kd, (k, 1, N)).doit()
    s_tbl = sum(table.values())
    s_orig = Sum(original, (k, lo, hi)).doit()

    print(f"\nSum over packed range (Piecewise)   = {s_pw}")
    print(f"Sum over packed range (Kronecker, symbolic Sum().doit()) = {s_kd}")
    print(f"Sum over table                      = {s_tbl}")
    print(f"Sum of original scattered summation = {s_orig}")

    assert s_pw == s_kd == s_tbl == s_orig, "sums disagree"
    for i in range(1, N + 1):
        assert pw.subs(k, i) == kd.subs(k, i) == lookup(i), f"mismatch at {i}"
    # packed values are exactly the stored ones in order, rest are 0
    assert [pw.subs(k, j + 1) for j in range(len(vals))] == vals
    assert all(pw.subs(k, i) == 0 for i in range(len(vals) + 1, N + 1))
    print("\nALL CHECKS PASSED: values moved to consecutive indices 1.."
          f"{len(vals)}, remaining {N - len(vals)} indices are 0, "
          "sum preserved, term stays symbolic (Sum + .subs both work).")
