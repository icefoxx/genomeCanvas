"""Global vertical layout allocation."""


def allocate_heights(preferred: list[float], total: float) -> list[float]:
    """Scale preferred track heights to a fixed total height."""
    if total <= 0 or any(height <= 0 for height in preferred):
        raise ValueError("heights and total must be positive")
    factor = total / sum(preferred)
    return [height * factor for height in preferred]
