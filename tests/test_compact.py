"""Pytest suite for the SymPy compact-summation prototype."""
import sympy as sp
from sympy import Sum, KroneckerDelta, symbols, S
from main import (stored_values, compact_piecewise,
                  compact_kronecker, compact_indexed)

k = symbols('k', integer=True)

# scattered term: non-zero only at 3, 17, 42, 99, 500
TERM = (7*KroneckerDelta(k, 3) + 4*KroneckerDelta(k, 17)
        + 9*KroneckerDelta(k, 42) - 2*KroneckerDelta(k, 99)
        + 5*KroneckerDelta(k, 500))
LO, HI = 1, 1000
VALS = [7, 4, 9, -2, 5]


def test_extraction_finds_scattered_values():
    got = stored_values(TERM, k, LO, HI)
    assert [i for i, _ in got] == [3, 17, 42, 99, 500]
    assert [v for _, v in got] == VALS


def test_values_packed_consecutively_from_1():
    pw = compact_piecewise(VALS, k, start=1)
    assert [pw.subs(k, j + 1) for j in range(len(VALS))] == VALS


def test_all_other_indexes_are_zero():
    pw = compact_piecewise(VALS, k, start=1)
    assert all(pw.subs(k, i) == 0 for i in range(len(VALS) + 1, HI + 1))


def test_three_forms_agree_pointwise():
    pw = compact_piecewise(VALS, k, start=1)
    kd = compact_kronecker(VALS, k, start=1)
    _, _, lookup = compact_indexed(VALS, start=1, length=HI)
    for i in range(1, HI + 1):
        assert pw.subs(k, i) == kd.subs(k, i) == lookup(i)


def test_kronecker_form_sums_symbolically():
    kd = compact_kronecker(VALS, k, start=1)
    assert Sum(kd, (k, 1, HI)).doit() == 23


def test_sum_preserved_across_compaction():
    before = Sum(TERM, (k, LO, HI)).doit()
    after = Sum(compact_kronecker(VALS, k, start=1), (k, 1, HI)).doit()
    assert before == after == 23


def test_configurable_start_offset():
    pw = compact_piecewise(VALS, k, start=10)
    assert [pw.subs(k, 10 + j) for j in range(len(VALS))] == VALS
    assert pw.subs(k, 9) == 0 and pw.subs(k, 15) == 0


def test_empty_input_is_all_zero():
    pw = compact_piecewise([], k, start=1)
    assert all(pw.subs(k, i) == 0 for i in range(1, 20))
