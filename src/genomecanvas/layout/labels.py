"""Priority-based label collision filtering."""


def visible_labels(
    labels: list[tuple[float, float, str, int]],
) -> list[tuple[float, float, str, int]]:
    """Keep higher-priority non-overlapping estimated text boxes."""
    accepted: list[tuple[float, float, str, int]] = []
    for label in sorted(labels, key=lambda value: (-value[3], value[0], value[1], value[2])):
        left, right, _, _ = label
        if all(right <= other[0] or left >= other[1] for other in accepted):
            accepted.append(label)
    return accepted


def label_layout(
    labels: list[tuple[float, float, str, int]], *, leader_lines: bool = False
) -> list[dict[str, object]]:
    """Return visible labels plus optional leader-line instructions."""
    visible = visible_labels(labels)
    return [
        {"left": left, "right": right, "text": text, "priority": priority,
         "leader_line": leader_lines}
        for left, right, text, priority in visible
    ]
