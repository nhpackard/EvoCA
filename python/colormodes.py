"""Display color modes: the single list shared by controls.py, display.py
and sdl_worker.py. Index = the `colormode` int passed to evoca_colorize()
and stored in .evoca recipes, so only ever APPEND (never reorder).
tests/test_colormodes.py asserts this list matches the C side's count."""

COLOR_MODES = [
    "genome",        # 0: LUT-hash color; full brightness v=1, dimmed v=0; dead black
    "env-food",      # 1
    "priv-food",     # 2
    "births",        # 3
    "age",           # 4
    "neighborhood",  # 5: 3x3 Moore pattern of v (9 bits, 512 colors; all-zero black)
    "food-waste",    # 6: green = env food F, red = env waste W (chem n_chem=2), blue = alive
]
