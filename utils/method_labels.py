"""Single source of truth for method identity across ALL figures and tables.

Internal method KEYS (e.g. "ESM", "DDG") are hard-coded in the survey and many
downstream scripts, so they MUST NOT change. What humans see, however, should be
the proper display name. Every plot/table that renders a method name should pull
it from here -- there is exactly ONE place to edit a label.

Usage:
    from utils.method_labels import METHOD_ORDER, disp, order_methods
    ax.set_xticklabels([disp(m) for m in cols])
"""

# Canonical left-to-right / top-to-bottom order for every cross-method plot.
METHOD_ORDER = ["AF2", "AF3", "ESM", "Boltz2", "DDG", "MSAT", "CCMpred", "S4PRED"]

# KEY -> human display string. Unicode (Δ = U+0394) renders in both matplotlib
# (DejaVu Sans) and HTML, so a single string serves every output target.
METHOD_DISPLAY = {
    "AF2":     "AF2",
    "AF3":     "AF3",
    "ESM":     "ESMFold2",
    "Boltz2":  "Boltz-2",
    "DDG":     "ΔΔG",      # ΔΔG
    "MSAT":    "MSA-Transformer",
    "CCMpred": "CCMpred",
    "S4PRED":  "S4PRED",
}


def disp(method: str) -> str:
    """Display label for a method key (falls back to the key itself)."""
    return METHOD_DISPLAY.get(method, method)


def order_methods(present):
    """Return `present` in canonical METHOD_ORDER, with any extras appended."""
    present = list(present)
    extra = [m for m in present if m not in METHOD_ORDER]
    return [m for m in METHOD_ORDER if m in present] + sorted(extra)
