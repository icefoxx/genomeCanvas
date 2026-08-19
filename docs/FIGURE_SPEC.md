# FigureSpec

FigureSpec supports explicit layout and scale controls:

```yaml
scales:
  mode: auto_highest       # default: shared maximum across signal tracks
  # mode: per_track        # independent automatic limits
  # mode: fixed
  # ymin: 0
  # ymax: 400
layout:
  width: 7.0
  height: 10.0
  aspect_ratio: 0.7       # width / height; use height or aspect_ratio
  track_gap: 0.0          # remove gaps between aligned tracks
  upper_annotation_ratio: 2.0  # total signal height : gene/cytoband height
  x_label: "GRCm38 chromosome 7 coordinate"
  show_x_axis: true
  show_gene_track: true
  show_cytoband: true
  cytoband_marker: bar     # bar, line, or triangle
  x_tick_count: 5
  y_axis_mode: all         # outer, all, left, right
  y_value_transform: identity  # identity, divide, log2, log10
  y_divisor: 10
  y_tick_count: 4
  y_tick_decimals: 1       # 1.42123 -> 1.4; 0 -> 1
  y_tick_label_pad: 3.0    # negative values move labels inside the plot
  show_y_ticks: true
  show_y_axis: true
  y_tick_length: 3.0
  show_major_frame: true
tracks:
  - id: input_rep1
    type: signal
    scale_group: input
    style:
      color: "#377eb8"
      ymin: 0
      ymax: 400
```

`auto_highest` uses the highest finite value in the shared scale group and adds
a small headroom margin. `fixed` values and per-track `style.ymin/ymax` override
automatic selection. No biological normalization is inferred.

Canvas shape is independently configurable. For a 16:10 canvas:

```yaml
layout:
  width: 16
  height: 10
  track_gap: 0.0
```
