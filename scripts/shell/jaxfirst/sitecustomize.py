# Auto-imported by Python at interpreter startup (it lives on PYTHONPATH and
# `site` imports any `sitecustomize` it finds). Purpose: force jaxlib to
# initialize its XLA backend and register xla_data.proto BEFORE TensorFlow is
# imported.
#
# Why: jaxlib and tensorflow each statically link a different XLA build and both
# register `xla/xla_data.proto` into protobuf's global C++ descriptor pool.
# jaxlib aborts ("File already exists in database: xla/xla_data.proto") when it
# is the SECOND registrant; TensorFlow tolerates being second. AlphaFold's
# embedding path (used by bioemu for the per-cluster MSA embedding) imports TF
# first and runs jax later -> jaxlib is second -> abort. Initializing jax here,
# at process start, makes jaxlib the FIRST registrant so TF's later import is
# the (tolerated) duplicate.
#
# Wrapped in try/except so a missing/broken jax never breaks unrelated python
# invocations that happen to have this dir on PYTHONPATH.
try:
    import jax.numpy as _jnp

    _jnp.zeros(1).block_until_ready()  # forces real XLA backend init, not just import
except Exception:
    pass
