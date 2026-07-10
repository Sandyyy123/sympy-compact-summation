# sympy-compact-summation

Move the **stored (non-zero) values** of a SymPy summation out of their scattered
indexes into a **consecutive block** at the start of a target range (e.g. `1..1000`),
packed next to each other until the values run out - every remaining index reads `0`.
The result stays a real SymPy object: `Sum()` and `.subs()` both keep working.

## The operation

| Stored value | Original (scattered) index | | New (packed) index |
|---|---|---|---|
| 7  | 3   | → | 1 |
| 4  | 17  | → | 2 |
| 9  | 42  | → | 3 |
| -2 | 99  | → | 4 |
| 5  | 500 | → | 5 |
| 0  | everything else | → | 6 .. 1000 |

Values are packed in ascending-index order and the summation total is preserved
(`7+4+9-2+5 = 23` before and after).

## Run it

```bash
pip install -r requirements.txt
python main.py
```

Expected tail of the output:

```
ALL CHECKS PASSED: values moved to consecutive indices 1..5, remaining 995
indices are 0, sum preserved, term stays symbolic (Sum + .subs both work).
```

## Three interchangeable representations

- **Piecewise** - value `j` at index `j`, else `0`. Matches the "indexes 1..1000" model directly.
- **KroneckerDelta** - pure algebra, no branching. This is the form that **sums
  symbolically** under `Sum(...).doit()`.
- **IndexedBase table** - an explicit `a[1..1000]` lookup; unset indexes read as `0`.

## A SymPy detail handled correctly

`Sum(Piecewise-with-Eq)` does **not** auto-collapse in SymPy, so the KroneckerDelta
form is used wherever the symbolic sum must actually evaluate. All three forms agree
index-by-index.

---
Prototype by Dr. Sandeep Grover.

## Tests

```bash
pip install pytest
python -m pytest -q tests/
```

8 tests cover extraction, consecutive packing, zero-fill, three-form agreement, symbolic sum, sum-preservation, configurable start, and empty input.
