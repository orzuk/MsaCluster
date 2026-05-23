"""SOV (Segment OVerlap) metric for secondary-structure similarity.

Pure-Python implementation of SOV'99 (Zemla et al. 1999, Proteins 34:220-223).
No external dependencies.

Quick reference
---------------
- SOV is the standard segment-level metric for evaluating SS predictions.
- Unlike Q3 (per-residue accuracy), SOV measures whether predicted SEGMENTS
  of H/E/C align with observed segments at the right positions, penalising
  broken-up predictions even when per-residue accuracy is high.
- Range 0..100 (or 0..1 if normalize=True). Higher = more similar.

For our pipeline
----------------
We compute SOV(predicted_ss_for_cluster, true_ss_for_fold) for each cluster
against BOTH folds of the pair. Signed difference SOV_F1 - SOV_F2 is the
per-cluster "SS-based" preference toward fold 1 vs fold 2.

Public API
----------
- segments(ss, states="HEC") -> list of (state, start, end) inclusive
- sov(ref, pred, states="HEC", per_state=False) -> overall SOV in 0..100
- sov_per_state(ref, pred) -> {state: SOV_state}

Both ref and pred are H/E/C strings of the SAME length. Characters not in
`states` are treated as "C" (coil) by default to match common practice.
"""

from __future__ import annotations


def _normalize_ss(s: str, states: str) -> str:
    """Map any character not in `states` to 'C' (coil) by default."""
    valid = set(states)
    return "".join(c if c in valid else "C" for c in s)


def segments(ss: str, states: str = "HEC") -> list[tuple[str, int, int]]:
    """Return list of (state, start_idx, end_idx) inclusive segments.
    Only segments whose state is in `states` are returned."""
    if not ss:
        return []
    ss = _normalize_ss(ss, states)
    out: list[tuple[str, int, int]] = []
    cur = ss[0]; start = 0
    for i in range(1, len(ss)):
        if ss[i] != cur:
            if cur in states:
                out.append((cur, start, i - 1))
            cur = ss[i]; start = i
    if cur in states:
        out.append((cur, start, len(ss) - 1))
    return out


def _overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Return length of overlap between two intervals [a0,a1] and [b0,b1].
    0 if no overlap."""
    lo = max(a[0], b[0])
    hi = min(a[1], b[1])
    return max(0, hi - lo + 1)


def sov_per_state(ref: str, pred: str, states: str = "HEC") -> dict[str, float]:
    """Return per-state SOV in 0..100. Returns NaN for states absent in ref."""
    if len(ref) != len(pred):
        raise ValueError(f"Length mismatch: ref={len(ref)} vs pred={len(pred)}")
    ref = _normalize_ss(ref, states)
    pred = _normalize_ss(pred, states)
    ref_segs = segments(ref, states)
    pred_segs = segments(pred, states)
    result: dict[str, float] = {}
    for state in states:
        ref_s = [s for s in ref_segs if s[0] == state]
        pred_s = [s for s in pred_segs if s[0] == state]
        N_state = sum(s[2] - s[1] + 1 for s in ref_s)
        if N_state == 0:
            result[state] = float("nan")
            continue
        total = 0.0
        for s1 in ref_s:
            len1 = s1[2] - s1[1] + 1
            # Find ALL predicted segments of this state that overlap s1
            overlapping = [s2 for s2 in pred_s
                           if _overlap((s1[1], s1[2]), (s2[1], s2[2])) > 0]
            if not overlapping:
                # No overlap: contribution is 0
                continue
            # Use ALL overlapping pairs (SOV'99 averaging convention).
            for s2 in overlapping:
                len2 = s2[2] - s2[1] + 1
                minov = _overlap((s1[1], s1[2]), (s2[1], s2[2]))
                maxov = max(s1[2], s2[2]) - min(s1[1], s2[1]) + 1
                delta = min(maxov - minov, minov, len1 // 2, len2 // 2)
                total += ((minov + delta) / maxov) * len1
        result[state] = 100.0 * total / N_state
    return result


def sov(ref: str, pred: str, states: str = "HEC",
         per_state: bool = False) -> float | dict[str, float]:
    """Overall SOV (weighted average across states by N_state).
    Returns 0..100. If per_state=True, returns the dict."""
    per = sov_per_state(ref, pred, states)
    if per_state:
        return per
    if len(ref) != len(pred):
        raise ValueError(f"Length mismatch: ref={len(ref)} vs pred={len(pred)}")
    ref = _normalize_ss(ref, states)
    total_N = sum(1 for c in ref if c in states)
    if total_N == 0:
        return float("nan")
    num = 0.0
    for state, v in per.items():
        if v != v:  # NaN
            continue
        N_state = sum(1 for c in ref if c == state)
        num += v * N_state
    return num / total_N


if __name__ == "__main__":
    # Sanity tests
    # Perfect prediction
    assert abs(sov("HHHHEEEECCCCHHHH", "HHHHEEEECCCCHHHH") - 100.0) < 1e-6
    # All wrong
    out = sov("HHHHEEEECCCCHHHH", "CCCCHHHHEEEECCCC")
    assert 0 <= out < 5, out
    # Random vs structured cases
    print("Perfect:", sov("HHHHEEEECCCCHHHH", "HHHHEEEECCCCHHHH"))
    print("Shift by 1:", sov("HHHHHHHH", "CHHHHHHC"))
    print("Per-state:", sov_per_state("HHHHEEEECCCCHHHH", "HHEEEEHHCCCCHHHH"))
