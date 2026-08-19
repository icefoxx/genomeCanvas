"""Small backend-neutral drawing instructions."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Line:
    x: tuple[float, ...]
    y: tuple[float, ...]
    color: str = "black"
    width: float = 1.0


@dataclass(frozen=True, slots=True)
class Rectangle:
    x: float
    y: float
    width: float
    height: float
    color: str = "black"
    alpha: float = 1.0


@dataclass(frozen=True, slots=True)
class Text:
    x: float
    y: float
    value: str
    size: float = 8.0
