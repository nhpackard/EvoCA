"""Which model a params dict or recipe describes.

The one place the `model` key and its load-time default live: a recipe or
params dict without `model` is EvoCA (every recipe and scan row written
before ChemCA existed), so the default is applied when reading, never
written back into old files. Dispatch points (evoca_explore.build_sim,
evoca_nb.make_sim, evoca_py.import_run) call model_of() and split_model().
"""

MODEL_KEY = 'model'
DEFAULT_MODEL = 'evoca'
MODELS = ('evoca', 'chemca')


def model_of(d):
    """The model named by dict d (params or recipe); 'evoca' if absent."""
    m = d.get(MODEL_KEY, DEFAULT_MODEL)
    if m not in MODELS:
        raise ValueError(f"unknown model {m!r}; known: {MODELS}")
    return m


def split_model(params):
    """(model, params without the model key)."""
    return model_of(params), {k: v for k, v in params.items() if k != MODEL_KEY}
