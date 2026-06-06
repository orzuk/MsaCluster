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
# Ordered roughly strongest/most-informative -> weakest: the structure
# predictors first, then DDG and S4PRED, then the sparse, noisy contact-based
# methods (MSA-Transformer, CCMpred) last.
METHOD_ORDER = ["AF2", "AF3", "ESM", "Boltz2", "BioEmu", "DDG", "S4PRED", "MSAT", "CCMpred"]

# KEY -> human display string. Unicode (Δ = U+0394) renders in both matplotlib
# (DejaVu Sans) and HTML, so a single string serves every output target.
METHOD_DISPLAY = {
    "AF2":     "AF2",
    "AF3":     "AF3",
    "ESM":     "ESMFold2",
    "Boltz2":  "Boltz-2",
    "BioEmu":  "BioEmu",
    "DDG":     "ΔΔG",      # ΔΔG
    "MSAT":    "MSA-Transformer",
    "CCMpred": "CCMpred",
    "S4PRED":  "S4PRED",
}


# KEY -> the NATIVE per-fold metric each method reports (what the TM1/TM2 slots
# of the unified survey actually hold). Used for per-method axis labels.
METHOD_METRIC = {
    "AF2":     "TM-score",
    "AF3":     "TM-score",
    "ESM":     "TM-score",
    "Boltz2":  "TM-score",
    "BioEmu":  "TM-score",
    "DDG":     "−ΣΔΔG",        # kcal/mol, higher = more stable on that fold
    "S4PRED":  "SSS",          # secondary-structure (H/E) composition similarity
    "MSAT":    "recall",       # contact-map recall
    "CCMpred": "recall",
}


def disp(method: str) -> str:
    """Display label for a method key (falls back to the key itself)."""
    return METHOD_DISPLAY.get(method, method)


def metric_label(method: str) -> str:
    """Native per-fold metric name for a method (for axis labels)."""
    return METHOD_METRIC.get(method, "support")


def order_methods(present):
    """Return `present` in canonical METHOD_ORDER, with any extras appended."""
    present = list(present)
    extra = [m for m in present if m not in METHOD_ORDER]
    return [m for m in METHOD_ORDER if m in present] + sorted(extra)
