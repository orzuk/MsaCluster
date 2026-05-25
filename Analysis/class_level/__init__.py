"""Class-level analysis: do trigger classes (mutation / ligand / oligomerization /
protein_binding / equilibrium_or_unknown) differ in per-cluster fold-preference
variability, aggregated across the 8 prediction methods?

The pair-level BH FDR analysis is severely under-powered (0 hits on full pool,
2 hits on stratified within-class BH at 85 pairs). This module asks a finer
question: even without per-pair FDR-significant hits, do trigger classes show
systematically different *rates* of clade-level signal?

Three complementary tests are implemented:
  - Fisher / Stouffer combination of 8 per-method Kruskal-Wallis p-values
    (combined_tests.py)
  - Linear mixed model with pair as random effect, method + class as fixed
    (lmm.py)
  - Permutation test on a multi-method aggregate F-ratio statistic
    (permutation.py)

Entry point: Analysis/class_level/run_all.py (or scripts/run_class_level_analysis.py).
"""

from .build_stats import build_per_pair_per_method_stats
from .combined_tests import per_method_kruskal, combine_methods
from .lmm import fit_class_lmm
from .permutation import permutation_test_class_effect

__all__ = [
    "build_per_pair_per_method_stats",
    "per_method_kruskal",
    "combine_methods",
    "fit_class_lmm",
    "permutation_test_class_effect",
]
