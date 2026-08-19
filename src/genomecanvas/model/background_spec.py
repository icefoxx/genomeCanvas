"""Declarative local-background transformation specifications."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, root_validator, validator


class BackgroundSpec(BaseModel):
    """Describe how local genomic background is estimated."""

    strategy: Literal["flanking", "rolling_window"] = "flanking"
    upstream: int | None = Field(default=None, ge=0)
    downstream: int | None = Field(default=None, ge=0)
    rolling_window: int | None = Field(default=None, ge=1)
    estimator: Literal["mean", "median", "trimmed_mean"] = "median"
    trim_fraction: float = Field(default=0.1, ge=0, lt=0.5)
    exclusions: list[dict[str, object]] = Field(default_factory=list)
    warning_threshold: float | None = Field(default=None, ge=0)

    @validator("upstream", "downstream")
    def nonzero_flank(cls, value: int) -> int:
        if value is not None and value < 0:
            raise ValueError("flank lengths must be non-negative")
        return value

    @root_validator(skip_on_failure=True)
    def strategy_parameters(cls, values: dict[str, object]) -> dict[str, object]:
        strategy = values.get("strategy")
        if strategy == "flanking" and (
            values.get("upstream") is None or values.get("downstream") is None
        ):
            raise ValueError("flanking strategy requires upstream and downstream lengths")
        if strategy == "rolling_window" and values.get("rolling_window") is None:
            raise ValueError("rolling_window strategy requires rolling_window")
        return values


class TransformationSpec(BaseModel):
    """A derived visualization transformation in FigureSpec."""

    id: str
    operation: Literal[
        "difference", "subtract", "ratio", "log2_ratio", "local_enrichment",
        "local_zscore", "smooth", "robust_scale",
    ]
    input: str
    background: BackgroundSpec | None = None
    pseudocount: float | None = Field(default=None, ge=0)
    label: str | None = None
