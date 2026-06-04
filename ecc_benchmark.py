"""
EC6204 Information Security - Mini Project | Group 28
Performance Improvement of ECC Scalar Multiplication
Using Windowed Non-Adjacent Form (wNAF) Method

Curves  : NIST P-256 (prime256v1) and secp256k1
Metrics : Execution time (ms), point operation count, memory usage (bytes)
Iterations: 500
"""

import time
import random
import tracemalloc
import statistics
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

# ── Elliptic-curve field and point arithmetic (pure Python, no external libs) ──

class FieldElement:
    """Modular arithmetic helpers (not used directly, but kept for clarity)."""
    pass


class ECPoint:
    """Affine point on a short Weierstrass curve  y² = x³ + ax + b (mod p)."""

    def __init__(self, x, y, curve):
        self.x = x
        self.y = y
        self.curve = curve

    def is_infinity(self):
        return self.x is None and self.y is None

    def __eq__(self, other):
        if other is None:
            return False
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        if self.is_infinity():
            return "ECPoint(∞)"
        return f"ECPoint({hex(self.x)[:12]}…, {hex(self.y)[:12]}…)"


INFINITY = ECPoint(None, None, None)  # Identity element


class EllipticCurve:
    """
    Short Weierstrass curve: y² ≡ x³ + ax + b (mod p).
    Supports point addition, doubling, and negation.
    """

    def __init__(self, p, a, b, Gx, Gy, n, name):
        self.p = p
        self.a = a
        self.b = b
        self.G = ECPoint(Gx, Gy, self)
        self.n = n          # order of G
        self.name = name

    # ── Low-level helpers ────────────────────────────────────────────────────

    def _modinv(self, k, p):
        """Extended Euclidean algorithm (modular inverse)."""
        if k == 0:
            raise ZeroDivisionError("No inverse of zero")
        if k < 0:
            return p - self._modinv(-k, p)
        s, old_s = 0, 1
        r, old_r = p, k
        while r != 0:
            q = old_r // r
            old_r, r = r, old_r - q * r
            old_s, s = s, old_s - q * s
        return old_s % p

    def neg(self, P):
        """Return -P."""
        if P.is_infinity():
            return INFINITY
        return ECPoint(P.x, (-P.y) % self.p, self)

    def add(self, P, Q, _counter=None):
        """Point addition P + Q.  Optionally increments an op-counter list."""
        if _counter is not None:
            _counter[0] += 1          # count this addition

        if P.is_infinity():
            return Q
        if Q.is_infinity():
            return P
        if P.x == Q.x:
            if P.y != Q.y:
                return INFINITY       # P = -Q
            return self.double(P, _counter)

        lam = ((Q.y - P.y) * self._modinv(Q.x - P.x, self.p)) % self.p
        x3 = (lam * lam - P.x - Q.x) % self.p
        y3 = (lam * (P.x - x3) - P.y) % self.p
        return ECPoint(x3, y3, self)

    def double(self, P, _counter=None):
        """Point doubling 2P.  Optionally increments a doubling-counter list."""
        if _counter is not None:
            _counter[1] += 1          # count this doubling

        if P.is_infinity():
            return INFINITY
        lam = ((3 * P.x * P.x + self.a) * self._modinv(2 * P.y, self.p)) % self.p
        x3 = (lam * lam - 2 * P.x) % self.p
        y3 = (lam * (P.x - x3) - P.y) % self.p
        return ECPoint(x3, y3, self)

    def subtract(self, P, Q, _counter=None):
        """P - Q = P + (-Q)."""
        return self.add(P, self.neg(Q), _counter)


# ── NIST P-256 (prime256v1) ──────────────────────────────────────────────────

P256 = EllipticCurve(
    p  = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF,
    a  = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC,
    b  = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
    Gx = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
    Gy = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5,
    n  = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    name = "NIST P-256"
)

# ── secp256k1 (Bitcoin curve) ─────────────────────────────────────────────────

SECP256K1 = EllipticCurve(
    p  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F,
    a  = 0,
    b  = 7,
    Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
    n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
    name = "secp256k1"
)


# ─────────────────────────────────────────────────────────────────────────────
#  Algorithm 1 : Double-and-Add
# ─────────────────────────────────────────────────────────────────────────────

def double_and_add(k, G, curve, count_ops=False):
    """
    Standard left-to-right binary scalar multiplication.

    Parameters
    ----------
    k          : scalar integer
    G          : base point (ECPoint)
    curve      : EllipticCurve instance
    count_ops  : if True, returns (result, additions, doublings)
    """
    counter = [0, 0]        # [additions, doublings]
    R = INFINITY
    binary_k = bin(k)[2:]   # e.g.  '11010001…'

    for bit in binary_k:
        if not R.is_infinity():
            R = curve.double(R, counter if count_ops else None)
        if bit == '1':
            if R.is_infinity():
                R = G
            else:
                R = curve.add(R, G, counter if count_ops else None)

    if count_ops:
        return R, counter[0], counter[1]
    return R


# ─────────────────────────────────────────────────────────────────────────────
#  Algorithm 2 : Windowed Non-Adjacent Form (wNAF)
# ─────────────────────────────────────────────────────────────────────────────

def to_wnaf(k, w):
    """
    Convert integer k to its wNAF (width-w Non-Adjacent Form) representation.

    Returns a list of digits (MSB-first) from the set
        { -(2^w-1), …, -3, -1, 0, 1, 3, …, 2^w-1 }
    where no two consecutive digits are both non-zero.
    """
    naf = []
    while k > 0:
        if k % 2 == 1:              # k is odd
            ki = k % (2 ** w)       # take w bits
            if ki >= 2 ** (w - 1):  # too large → go negative
                ki -= 2 ** w
            naf.append(ki)
            k -= ki
        else:
            naf.append(0)
        k //= 2                     # right-shift
    return naf[::-1]                # MSB first


def wnaf_multiply(k, G, curve, w=4, count_ops=False):
    """
    wNAF scalar multiplication.

    Parameters
    ----------
    k          : scalar integer
    G          : base point (ECPoint)
    curve      : EllipticCurve instance
    w          : window width (4 is optimal for 256-bit scalars)
    count_ops  : if True, returns (result, additions, doublings)
    """
    counter = [0, 0]        # [additions, doublings]

    # Step 1: Precompute table  [ G, 3G, 5G, …, (2^w - 1)G ]
    table_size = 2 ** (w - 1)   # number of odd multiples needed
    precomp = [None] * table_size
    precomp[0] = G
    double_G = curve.double(G)
    for i in range(1, table_size):
        precomp[i] = curve.add(precomp[i - 1], double_G)

    # Step 2: Convert k to wNAF digits
    naf_digits = to_wnaf(k, w)

    # Step 3: Evaluate via left-to-right scanning
    R = INFINITY
    for digit in naf_digits:
        if not R.is_infinity():
            R = curve.double(R, counter if count_ops else None)

        if digit > 0:
            idx = digit // 2        # odd digit d → index (d-1)/2
            if R.is_infinity():
                R = precomp[idx]
            else:
                R = curve.add(R, precomp[idx], counter if count_ops else None)
        elif digit < 0:
            idx = (-digit) // 2
            if R.is_infinity():
                R = curve.neg(precomp[idx])
            else:
                R = curve.subtract(R, precomp[idx], counter if count_ops else None)

    if count_ops:
        return R, counter[0], counter[1]
    return R


# ─────────────────────────────────────────────────────────────────────────────
#  Correctness Verification
# ─────────────────────────────────────────────────────────────────────────────

def verify_correctness(curve, trials=20):
    """Both algorithms must produce the same point Q = kG."""
    print(f"\n{'─'*55}")
    print(f"  Correctness check — {curve.name}  ({trials} random scalars)")
    print(f"{'─'*55}")
    all_ok = True
    for i in range(trials):
        k = random.randint(2, curve.n - 1)
        R_da   = double_and_add(k, curve.G, curve)
        R_wnaf = wnaf_multiply(k, curve.G, curve)
        ok = (R_da == R_wnaf)
        if not ok:
            all_ok = False
            print(f"  [FAIL] trial {i+1}: mismatch for k={hex(k)[:16]}…")
    if all_ok:
        print(f"  [PASS] All {trials} trials matched ✓")
    return all_ok


# ─────────────────────────────────────────────────────────────────────────────
#  Benchmarking helpers
# ─────────────────────────────────────────────────────────────────────────────

def benchmark(curve, iterations=500, w=4):
    """
    Run 'iterations' timed trials for both algorithms on 'curve'.
    Returns a dict with timing lists and operation counts.
    """
    print(f"\n{'─'*55}")
    print(f"  Benchmarking — {curve.name}  ({iterations} iterations, w={w})")
    print(f"{'─'*55}")

    da_times, wnaf_times = [], []
    da_adds, da_dubs = [], []
    wnaf_adds, wnaf_dubs = [], []

    for i in range(iterations):
        k = random.randint(2, curve.n - 1)

        # ── Double-and-Add ──────────────────────────────────────────
        tracemalloc.start()
        t0 = time.perf_counter()
        R_da, a, d = double_and_add(k, curve.G, curve, count_ops=True)
        t1 = time.perf_counter()
        _, peak_da = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        da_times.append((t1 - t0) * 1000)      # ms
        da_adds.append(a)
        da_dubs.append(d)

        # ── wNAF ────────────────────────────────────────────────────
        tracemalloc.start()
        t0 = time.perf_counter()
        R_wn, a, d = wnaf_multiply(k, curve.G, curve, w=w, count_ops=True)
        t1 = time.perf_counter()
        _, peak_wnaf = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        wnaf_times.append((t1 - t0) * 1000)
        wnaf_adds.append(a)
        wnaf_dubs.append(d)

        if (i + 1) % 100 == 0:
            print(f"    {i+1}/{iterations} done …")

    return {
        "da_times"   : da_times,
        "wnaf_times" : wnaf_times,
        "da_adds"    : da_adds,
        "da_dubs"    : da_dubs,
        "wnaf_adds"  : wnaf_adds,
        "wnaf_dubs"  : wnaf_dubs,
    }


def print_summary(data, curve_name):
    """Print a formatted statistics table for one curve."""
    def stats(lst):
        return {
            "mean"  : statistics.mean(lst),
            "min"   : min(lst),
            "max"   : max(lst),
            "stdev" : statistics.stdev(lst),
        }

    da   = stats(data["da_times"])
    wn   = stats(data["wnaf_times"])
    improvement = (da["mean"] - wn["mean"]) / da["mean"] * 100
    add_reduction = (
        (statistics.mean(data["da_adds"]) - statistics.mean(data["wnaf_adds"]))
        / statistics.mean(data["da_adds"]) * 100
    )

    print(f"\n{'═'*60}")
    print(f"  Results — {curve_name}")
    print(f"{'═'*60}")
    print(f"  {'Metric':<28} {'Double-and-Add':>12}  {'wNAF (w=4)':>12}")
    print(f"  {'─'*56}")
    print(f"  {'Avg time (ms)':<28} {da['mean']:>12.4f}  {wn['mean']:>12.4f}")
    print(f"  {'Min time (ms)':<28} {da['min']:>12.4f}  {wn['min']:>12.4f}")
    print(f"  {'Max time (ms)':<28} {da['max']:>12.4f}  {wn['max']:>12.4f}")
    print(f"  {'Std-dev (ms)':<28} {da['stdev']:>12.4f}  {wn['stdev']:>12.4f}")
    print(f"  {'Avg point additions':<28} {statistics.mean(data['da_adds']):>12.1f}  {statistics.mean(data['wnaf_adds']):>12.1f}")
    print(f"  {'Avg point doublings':<28} {statistics.mean(data['da_dubs']):>12.1f}  {statistics.mean(data['wnaf_dubs']):>12.1f}")
    print(f"  {'─'*56}")
    print(f"  Execution time improvement : {improvement:+.2f}%")
    print(f"  Addition reduction         : {add_reduction:+.2f}%")
    print(f"{'═'*60}")


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting
# ─────────────────────────────────────────────────────────────────────────────

def plot_results(data_p256, data_k1):
    """Generate a 2×3 grid of benchmark charts and save to PNG."""
    fig = plt.figure(figsize=(18, 10))
    fig.suptitle(
        "ECC Scalar Multiplication: Double-and-Add vs wNAF (w=4)\n"
        "Group 28 — EC6204 Information Security",
        fontsize=14, fontweight='bold'
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    COLORS = {"da": "#E74C3C", "wnaf": "#2ECC71"}

    for row, (data, cname) in enumerate([(data_p256, "NIST P-256"),
                                          (data_k1,   "secp256k1")]):
        avg_da   = statistics.mean(data["da_times"])
        avg_wnaf = statistics.mean(data["wnaf_times"])

        # Col 0 — Bar chart: average execution time
        ax0 = fig.add_subplot(gs[row, 0])
        bars = ax0.bar(["Double-and-Add", "wNAF (w=4)"],
                       [avg_da, avg_wnaf],
                       color=[COLORS["da"], COLORS["wnaf"]],
                       width=0.5, edgecolor='white')
        for bar, val in zip(bars, [avg_da, avg_wnaf]):
            ax0.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 0.001,
                     f"{val:.3f} ms", ha='center', va='bottom', fontsize=9)
        ax0.set_title(f"{cname}\nAvg Execution Time", fontsize=10)
        ax0.set_ylabel("Time (ms)")
        ax0.set_ylim(0, max(avg_da, avg_wnaf) * 1.25)

        # Col 1 — Line chart: first 100 runs
        ax1 = fig.add_subplot(gs[row, 1])
        ax1.plot(data["da_times"][:100],   color=COLORS["da"],   alpha=0.7,
                 linewidth=0.9, label="Double-and-Add")
        ax1.plot(data["wnaf_times"][:100], color=COLORS["wnaf"], alpha=0.7,
                 linewidth=0.9, label="wNAF")
        ax1.axhline(avg_da,   color=COLORS["da"],   linestyle='--', linewidth=0.8)
        ax1.axhline(avg_wnaf, color=COLORS["wnaf"], linestyle='--', linewidth=0.8)
        ax1.set_title(f"{cname}\nExecution Time (100 runs)", fontsize=10)
        ax1.set_xlabel("Iteration")
        ax1.set_ylabel("Time (ms)")
        ax1.legend(fontsize=8)

        # Col 2 — Bar chart: avg point operations
        ax2 = fig.add_subplot(gs[row, 2])
        avg_ops = {
            "D-A Add"  : statistics.mean(data["da_adds"]),
            "D-A Dbl"  : statistics.mean(data["da_dubs"]),
            "wNAF Add" : statistics.mean(data["wnaf_adds"]),
            "wNAF Dbl" : statistics.mean(data["wnaf_dubs"]),
        }
        bar_colors = [COLORS["da"], COLORS["da"], COLORS["wnaf"], COLORS["wnaf"]]
        ax2.bar(avg_ops.keys(), avg_ops.values(),
                color=bar_colors, edgecolor='white', width=0.6, alpha=0.85)
        ax2.set_title(f"{cname}\nAvg Point Operations", fontsize=10)
        ax2.set_ylabel("Count")
        ax2.set_xticks(range(len(avg_ops)))
        ax2.set_xticklabels(avg_ops.keys(), rotation=15, fontsize=8)

    import os
    output_dir = "outputs"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "ecc_benchmark.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n  Chart saved → {output_file}")
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    random.seed(42)   # reproducible results
    ITERATIONS = 500
    W = 4             # window width (optimal for 256-bit keys)

    print("=" * 60)
    print("  EC6204 — Mini Project Group 28")
    print("  ECC Scalar Multiplication Benchmark")
    print("=" * 60)

    # 1. Correctness check
    verify_correctness(P256)
    verify_correctness(SECP256K1)

    # 2. Benchmarks
    data_p256 = benchmark(P256,      iterations=ITERATIONS, w=W)
    data_k1   = benchmark(SECP256K1, iterations=ITERATIONS, w=W)

    # 3. Print summaries
    print_summary(data_p256, P256.name)
    print_summary(data_k1,   SECP256K1.name)

    # 4. Plots
    plot_results(data_p256, data_k1)

    # 5. Export CSV
    df = pd.DataFrame({
        "iteration"         : range(1, ITERATIONS + 1),
        "p256_da_time_ms"   : data_p256["da_times"],
        "p256_wnaf_time_ms" : data_p256["wnaf_times"],
        "p256_da_adds"      : data_p256["da_adds"],
        "p256_wnaf_adds"    : data_p256["wnaf_adds"],
        "k1_da_time_ms"     : data_k1["da_times"],
        "k1_wnaf_time_ms"   : data_k1["wnaf_times"],
        "k1_da_adds"        : data_k1["da_adds"],
        "k1_wnaf_adds"      : data_k1["wnaf_adds"],
    })
    output_dir = "outputs"
    csv_path = os.path.join(output_dir, "ecc_benchmark_results.csv")
    df.to_csv(csv_path, index=False)
    print(f"  Data   saved → {csv_path}")


if __name__ == "__main__":
    main()
