"""Publication-oriented dimension helpers."""


PRESETS_INCH = {"single-column": 3.5, "1.5-column": 5.5, "double-column": 7.0}


def width_inches(name: str) -> float:
    """Resolve a named width preset."""
    if name not in PRESETS_INCH:
        raise ValueError(f"unknown width preset: {name}")
    return PRESETS_INCH[name]
