'use strict';
/* GenomeCanvas Studio front end: projects, data, declarative settings, jobs, saved figures. */

const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
function h(tag, attrs = {}, ...kids) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v == null || v === false) continue;
    if (k === 'class') el.className = v;
    else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else if (k === 'html') el.innerHTML = v;
    else if (k === 'value') el.value = v;
    else if (k === 'checked') el.checked = !!v;
    else el.setAttribute(k, v === true ? '' : v);
  }
  for (const kid of kids.flat()) {
    if (kid == null || kid === false) continue;
    el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  }
  return el;
}
/** append() stringifies null; skip absent optional children instead. */
const put = (el, ...kids) => { el.append(...kids.flat().filter(k => k != null && k !== false)); return el; };
/** replaceChildren() also stringifies null, so conditional children go through this. */
const fill = (el, ...kids) => { el.replaceChildren(...kids.flat().filter(k => k != null && k !== false)); return el; };
const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
const fmtBytes = n => n < 1024 ? n + ' B' : n < 1048576 ? (n / 1024).toFixed(0) + ' KB' : n < 1073741824 ? (n / 1048576).toFixed(1) + ' MB' : (n / 1073741824).toFixed(2) + ' GB';
const fmtBp = n => Number(n).toLocaleString('en-US');
const fmtDate = t => t ? new Date(t * 1000).toLocaleString(undefined, { dateStyle: 'medium', timeStyle: 'short' }) : '';
const slug = s => String(s || 'figure').toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '').slice(0, 60) || 'figure';
const normChrom = c => String(c || '').replace(/^chr/i, '').toLowerCase();

async function api(path, opt = {}) {
  const headers = { 'X-GenomeCanvas': 'studio', ...(typeof opt.body === 'string' ? { 'Content-Type': 'application/json' } : {}), ...opt.headers };
  const res = await fetch(path, { ...opt, headers, credentials: 'same-origin' });
  if (res.status === 401) { location.href = '/login'; throw new Error('Please sign in'); }
  let data;
  try { data = await res.json(); } catch { data = { error: res.statusText }; }
  if (!res.ok) throw new Error(data.error || 'Request failed');
  return data;
}
const papi = (sub = '', opt) => api(`/api/projects/${encodeURIComponent(S.project.id)}${sub}`, opt);

/* ------------------------------------------------------------------ state */
const STORE = 'gc-studio-v3';
const S = {
  meta: null, me: null, projects: [], project: null, role: 'editor', templates: [], figures: [],
  category: 1, plot: 2, files: [], settings: {}, advanced: '', live: false,
  zoom: 'fit', checker: false, history: [], current: null, revision: 0, resultRevision: -1,
  busy: false, queued: false, openFile: null, figureId: null, view: 'interactive',
  sections: { color: true, text: true, region: true, matrix: true, heatmap: true, profile: true, analysis: true },
  search: '',
};
const readOnly = () => S.role === 'viewer';
const store = {
  get(key, fallback) { try { return JSON.parse(localStorage.getItem(key)) ?? fallback; } catch { return fallback; } },
  set(key, value) { try { localStorage.setItem(key, JSON.stringify(value)); } catch { /* storage unavailable */ } },
};
const EDITABLE = ['role', 'label', 'group', 'color', 'cmap'];
/** Only user edits are remembered, so manifest updates by the core still flow through. */
function fileOverrides(f) {
  const out = { enabled: f.enabled !== false };
  for (const k of EDITABLE) if ((f[k] || '') !== (f.base[k] || '')) out[k] = f[k];
  return out;
}
function persist() {
  store.set(STORE, { live: S.live, sections: S.sections, project: S.project?.id });
  if (!S.project) return;
  store.set(`${STORE}:${S.project.id}`, {
    category: S.category, plot: S.plot, settings: S.settings, advanced: S.advanced,
    files: Object.fromEntries(S.files.map(f => [f.id, fileOverrides(f)])), order: S.files.map(f => f.id),
  });
}

function status(msg, kind) { const el = $('#status'); el.textContent = msg; el.className = kind === 'error' ? 'err' : ''; }
let toastTimer;
function toast(msg) { const t = $('#toast'); t.textContent = msg; t.hidden = false; clearTimeout(toastTimer); toastTimer = setTimeout(() => { t.hidden = true; }, 3600); }

/** Small modal form. fields: [{name,label,type:'text'|'textarea'|'select'|'check',value,options,help}] */
function ask({ title, fields, ok = 'Save', html }) {
  return new Promise(resolve => {
    const dlg = $('#dialog');
    $('#dialogtitle').textContent = title;
    $('#dialogok').textContent = ok;
    const inputs = {};
    fill($('#dialogbody'), ...(html ? [html] : []), ...(fields || []).map(f => {
      let input;
      if (f.type === 'textarea') input = h('textarea', { rows: 3, value: f.value || '' });
      else if (f.type === 'select') input = h('select', {}, ...f.options.map(([v, t]) => h('option', { value: v, selected: v === f.value }, t)));
      else if (f.type === 'check') input = h('input', { type: 'checkbox', checked: !!f.value });
      else input = h('input', { type: f.type || 'text', value: f.value ?? '', required: f.required, placeholder: f.placeholder || '' });
      inputs[f.name] = input;
      return f.type === 'check'
        ? h('label', { class: 'dialog-check' }, input, f.label)
        : h('label', { class: 'dialog-field' }, h('span', {}, f.label), input, f.help ? h('small', {}, f.help) : null);
    }));
    $('#dialogok').hidden = ok === null;
    dlg.onclose = () => {
      if (dlg.returnValue !== 'ok') { resolve(null); return; }
      resolve(Object.fromEntries(Object.entries(inputs).map(([k, el]) => [k, el.type === 'checkbox' ? el.checked : el.value.trim()])));
    };
    dlg.returnValue = '';
    dlg.showModal();
    Object.values(inputs)[0]?.focus();
  });
}

/* ------------------------------------------------------------------ catalog helpers */
const plotName = (c = S.category, p = S.plot) => S.meta.catalog[c][1][p];
const enabled = () => S.files.filter(f => f.enabled !== false);
const byRole = role => enabled().filter(f => f.role === role);
function paletteColors() {
  const custom = String(val('custom_palette') || '').split(/[,;\s]+/).filter(c => /^#?[0-9a-f]{3,8}$/i.test(c)).map(c => c.startsWith('#') ? c : '#' + c);
  if (custom.length) return custom;
  return (S.meta.palettes[val('palette')] || S.meta.palettes.studio).colors;
}
/** Colour the server will use: explicit per-file colour, else palette slot (signals first, then regions). */
function fileColor(f) {
  if (f.color) return f.color;
  if (f.enabled === false || f.role === 'genes') return null;
  const colors = paletteColors();
  const sig = byRole('signal');
  const idx = f.role === 'signal' ? sig.indexOf(f) : sig.length + byRole('regions').indexOf(f);
  return colors[idx % colors.length];
}
function requirement(c = S.category, p = S.plot) {
  const need = { ...S.meta.requirements[c][p] };
  if (c === 2 && val('rel_mode') !== 'regions') delete need.regions;
  const have = { signal: byRole('signal').length, regions: byRole('regions').length, genes: byRole('genes').length };
  const items = Object.entries(need).map(([role, n]) => ({ role, n, have: have[role], ok: have[role] >= n }));
  return { items, ok: items.every(i => i.ok) };
}
const ROLE_NAMES = { signal: ['signal', 'signals'], regions: ['region set', 'region sets'], genes: ['gene annotation', 'gene annotations'] };

/* ------------------------------------------------------------------ settings schema */
const C = (...cats) => c => cats.includes(c);
const matrixWhen = (c, p) => c === 1 || (c === 3 && p === 1);
const genomeWide = c => c === 2 && val('rel_mode') !== 'regions';
const regionWhen = c => c === 0 || c === 4 || (genomeWide(c) && val('rel_scope') === 'locus');
const DEFAULT_SKIP = 'chrM, MT, *_random, chrUn*, *_alt, *_fix, chrEBV, *hap*, *decoy*';
const heatWhen = (c, p) => c === 1 && p !== 1;
const profWhen = (c, p) => c === 1 && (p === 1 || p === 2);
const ASSEMBLIES = ['hg38', 'hg19', 'mm39', 'mm10', 'rn7', 'dm6', 'dm3', 'ce11', 'ce10', 'danRer11', 'sacCer3', 'TAIR10'];
const LEGEND = ['best', 'upper right', 'upper left', 'lower right', 'lower left', 'center right', 'center left', 'upper center', 'lower center', 'center'];
const unit = () => (val('size_unit') === 'mm' ? 'mm' : 'in');
const SECTIONS = [
  ['region', 'Region'], ['matrix', 'Regions & binning'], ['heatmap', 'Heatmap'], ['profile', 'Profile'],
  ['analysis', 'Analysis'], ['tracks', 'Track layout'], ['color', 'Colours & style'], ['text', 'Text & labels'],
  ['figure', 'Figure & export'], ['advanced', 'Advanced'],
];
const FIELDS = [
  // Region
  { k: 'chrom', l: 'Chromosome', t: 'combo', sec: 'region', d: 'chr1', list: () => chromList(), when: regionWhen },
  { k: 'assembly', l: 'Assembly', t: 'combo', sec: 'region', d: () => S.project?.genome || '', ph: 'e.g. hg38', list: ASSEMBLIES, when: C(0) },
  { k: 'start', l: 'Start (0-based)', t: 'number', sec: 'region', d: 0, min: 0, step: 1, when: regionWhen },
  { k: 'end', l: 'End (exclusive)', t: 'number', sec: 'region', d: 100000, min: 1, step: 1, when: regionWhen },
  { k: '_locus', t: 'locus', sec: 'region', wide: 1, when: regionWhen },
  { k: 'highlights', l: 'Highlight windows', t: 'textarea', sec: 'region', wide: 1, d: '', ph: '31136000-31138000 #f5c542 0.25', help: 'One window per line: start-end [colour] [opacity]. Coordinates may include commas or a chr: prefix.', when: C(0) },
  // Regions & binning
  { k: 'alignment', l: 'Anchor', t: 'select', sec: 'matrix', d: 'center', o: [['center', 'Region center'], ['TSS', 'Start / TSS'], ['TES', 'End / TES'], ['scale-regions', 'Scale regions']], when: matrixWhen },
  { k: 'flank', l: 'Flank (bp)', t: 'number', sec: 'matrix', d: 1000, min: 20, step: 50, when: matrixWhen },
  { k: 'binsize', l: 'Bin size (bp)', t: 'number', sec: 'matrix', d: 50, min: 1, step: 5, when: matrixWhen },
  { k: 'body_length', l: 'Scaled body (bp)', t: 'number', sec: 'matrix', d: 1000, min: 10, step: 100, when: (c, p) => matrixWhen(c, p) && val('alignment') === 'scale-regions' },
  { k: 'skip_zeros', l: 'Skip all-zero regions', t: 'check', sec: 'matrix', d: true, when: matrixWhen },
  // Heatmap
  { k: 'sort', l: 'Row order', t: 'select', sec: 'heatmap', d: 'descend', o: [['descend', 'Descending'], ['ascend', 'Ascending'], ['none', 'As in file']], when: heatWhen },
  { k: 'sort_by', l: 'Sort by', t: 'select', sec: 'heatmap', d: 'mean', o: ['mean', 'median', 'max', 'min', 'sum'], when: (c, p) => heatWhen(c, p) && val('sort') !== 'none' },
  { k: 'cluster', l: 'Clustering', t: 'select', sec: 'heatmap', d: (c, p) => (p === 3 ? 'kmeans' : 'none'), o: [['none', 'None'], ['kmeans', 'k-means'], ['hclust', 'Hierarchical']], when: heatWhen },
  { k: 'clusters', l: 'Clusters (k)', t: 'number', sec: 'heatmap', d: 3, min: 2, max: 20, step: 1, when: (c, p) => heatWhen(c, p) && val('cluster') !== 'none' },
  { k: 'normalization', l: 'Normalisation', t: 'select', sec: 'heatmap', d: 'raw', o: [['raw', 'None (raw)'], ['profileMinMax', 'Min–max per sample'], ['qnorm', 'Quantile']], when: heatWhen },
  { k: 'autoscale', l: 'Automatic colour range', t: 'check', sec: 'heatmap', d: true, when: heatWhen },
  { k: 'zmin', l: 'Colour min', t: 'number', sec: 'heatmap', d: '', ph: 'auto', step: 'any', when: heatWhen, dim: () => val('autoscale') },
  { k: 'zmax', l: 'Colour max', t: 'number', sec: 'heatmap', d: '', ph: 'auto', step: 'any', when: heatWhen, dim: () => val('autoscale') },
  { k: 'indi', l: 'Separate colour key per sample', t: 'check', sec: 'heatmap', d: true, when: heatWhen },
  { k: 'show_colorbar', l: 'Show colour key', t: 'check', sec: 'heatmap', d: true, when: heatWhen },
  { k: 'show_group_labels', l: 'Show group strip', t: 'check', sec: 'heatmap', d: true, when: heatWhen },
  { k: 'interpolation', l: 'Pixel smoothing', t: 'combo', sec: 'heatmap', d: 'bilinear', list: ['bilinear', 'nearest', 'bicubic', 'gaussian', 'spline16', 'none'], when: heatWhen },
  { k: 'width_per_sample', l: 'Panel width (in)', t: 'range', sec: 'heatmap', d: 2.7, min: 1, max: 6, step: 0.1, when: heatWhen },
  // Profile
  { k: 'plot_scale', l: 'Signal scale', t: 'select', sec: 'profile', d: 'linear', o: [['linear', 'Linear'], ['log2', 'log₂(x + 1)'], ['zscore', 'z-score']], when: C(1) },
  { k: 'x_unit', l: 'Distance unit', t: 'select', sec: 'profile', d: 'auto', o: ['auto', 'bp', 'kb', 'Mb'], when: C(1) },
  { k: 'error_estimates', l: 'Show SEM ribbon', t: 'check', sec: 'profile', d: true, when: profWhen },
  { k: 'ln_v', l: 'Vertical guide at anchor', t: 'check', sec: 'profile', d: true, when: C(1) },
  { k: 'ln_h', l: 'Horizontal guide line', t: 'check', sec: 'profile', d: false, when: (c, p) => c === 1 && p === 1 },
  { k: 'ln_h_pos', l: 'Guide at Y =', t: 'number', sec: 'profile', d: 0, step: 'any', when: (c, p) => c === 1 && p === 1 && val('ln_h') },
  { k: 'ylim_manual', l: 'Fixed Y range', t: 'check', sec: 'profile', d: false, when: profWhen },
  { k: 'ymin', l: 'Y min', t: 'number', sec: 'profile', d: '', step: 'any', when: (c, p) => profWhen(c, p) && val('ylim_manual') },
  { k: 'ymax', l: 'Y max', t: 'number', sec: 'profile', d: '', step: 'any', when: (c, p) => profWhen(c, p) && val('ylim_manual') },
  { k: 'legend', l: 'Show legend', t: 'check', sec: 'profile', d: true, when: (c, p) => c === 1 && p === 1 },
  { k: 'legend_pos', l: 'Legend position', t: 'select', sec: 'profile', d: 'best', o: LEGEND, when: (c, p) => c === 1 && p === 1 && val('legend') },
  // Analysis: what each sample is summarised over (deepTools multiBigwigSummary BED-file / bins)
  { k: 'rel_mode', l: 'Summarise signal over', t: 'select', sec: 'analysis', d: 'regions', wide: 1, o: [['regions', 'Selected regions (BED files)'], ['bins', 'Genome-wide bins'], ['windows', 'Sliding windows']], help: 'Genome-wide modes ignore region files and compare samples across the whole genome, like multiBigwigSummary bins.', when: C(2) },
  { k: 'rel_bin', l: 'Bin / window size (bp)', t: 'combo', sec: 'analysis', d: 10000, list: ['500', '1000', '2000', '5000', '10000', '50000', '100000'], when: genomeWide },
  { k: 'rel_step', l: 'Step (bp)', t: 'combo', sec: 'analysis', d: 2000, list: ['500', '1000', '2000', '5000'], help: 'Window size must be a multiple of the step.', when: c => genomeWide(c) && val('rel_mode') === 'windows' },
  { k: 'rel_scope', l: 'Scope', t: 'select', sec: 'analysis', d: 'genome', o: [['genome', 'Whole genome'], ['chroms', 'Selected chromosomes'], ['locus', 'Current region only']], when: genomeWide },
  { k: 'rel_chroms', l: 'Chromosomes', t: 'text', sec: 'analysis', d: '', ph: 'chr1, chr2, chrX', wide: 1, when: c => genomeWide(c) && val('rel_scope') === 'chroms' },
  { k: 'rel_skip', l: 'Skip chromosomes (patterns)', t: 'text', sec: 'analysis', d: DEFAULT_SKIP, wide: 1, help: 'Comma-separated, * wildcards; matched case-insensitively.', when: c => genomeWide(c) && val('rel_scope') !== 'locus' },
  { k: 'rel_blacklist', l: 'Blacklist', t: 'select', sec: 'analysis', d: '', wide: 1, o: () => [['', 'None'], ...S.files.filter(f => f.role === 'regions').map(f => [f.id, f.label])], help: 'Bins overlapping these regions are removed (e.g. ENCODE blacklist). Any region file in the project can be used.', when: genomeWide },
  { k: 'rel_stat', l: 'Statistic', t: 'select', sec: 'analysis', d: () => (val('rel_mode') === 'regions' ? 'sum' : 'mean'), o: [['mean', 'Mean (covered bases)'], ['mean0', 'Mean (missing = 0)'], ['sum', 'Sum (signal × bp)'], ['max', 'Maximum']], when: C(2) },
  { k: 'rel_skip_zeros', l: 'Skip all-zero rows', t: 'check', sec: 'analysis', d: true, when: C(2) },
  { k: 'rel_max_rows', l: 'Subsample to (rows)', t: 'number', sec: 'analysis', d: '', ph: 'all', min: 100, step: 1000, help: 'Random, reproducible subset for speed (seed 7201).', when: C(2) },
  { k: 'corr_transform', l: 'Transform before correlating', t: 'select', sec: 'analysis', d: 'raw', o: [['raw', 'None'], ['log1p', 'log(x + 1)']], help: 'log(x+1) or Spearman reduce the influence of a few very high bins.', when: (c, p) => c === 2 && p === 0 },
  { k: 'pca_top', l: 'Most variable features', t: 'number', sec: 'analysis', d: 1000, min: 2, step: 500, when: (c, p) => c === 2 && p === 1 },
  { k: 'method', l: 'Correlation', t: 'select', sec: 'analysis', d: 'pearson', o: [['pearson', 'Pearson'], ['spearman', 'Spearman']], when: (c, p) => c === 2 && p === 0 },
  { k: 'annotate', l: 'Print coefficients in cells', t: 'check', sec: 'analysis', d: true, when: (c, p) => c === 2 && p === 0 },
  { k: 'log', l: 'Signal transform', t: 'select', sec: 'analysis', d: 'log2', o: [['log2', 'log₂'], ['raw', 'Raw']], when: (c, p) => c === 2 && p === 1 },
  { k: 'point_size', l: 'Point size', t: 'range', sec: 'analysis', d: 60, min: 10, max: 300, step: 5, when: (c, p) => c === 2 && p === 1 },
  { k: 'maxrows', l: 'Max sampled regions', t: 'number', sec: 'analysis', d: 1000, min: 10, max: 20000, step: 100, when: (c, p) => c === 2 && p === 2 },
  { k: 'numerator', l: 'Numerator group', t: 'combo', sec: 'analysis', d: 'Treatment', list: () => groupList(), when: (c, p) => c === 3 && p === 0 },
  { k: 'denominator', l: 'Denominator group', t: 'combo', sec: 'analysis', d: 'Control', list: () => groupList(), when: (c, p) => c === 3 && p === 0 },
  { k: '_groups', t: 'groups', sec: 'analysis', wide: 1, when: (c, p) => c === 3 && p === 0 },
  { k: 'fdr', l: 'FDR threshold', t: 'number', sec: 'analysis', d: 0.05, min: 0, max: 1, step: 0.01, when: (c, p) => c === 3 && p === 0 },
  { k: 'center_bins', l: 'Center half-window (bins)', t: 'number', sec: 'analysis', d: 5, min: 1, max: 200, step: 1, when: (c, p) => c === 3 && p === 1 },
  { k: 'direction', l: 'Signal order', t: 'select', sec: 'analysis', d: 'b-a', o: [['b-a', 'A = first, B = second'], ['a-b', 'Swap A and B']], help: 'A and B are the first two enabled signals. Drag files to reorder.', wide: 1, when: C(4) },
  { k: 'pseudocount', l: 'Pseudocount', t: 'number', sec: 'analysis', d: 0.01, min: 0.000001, step: 'any', when: (c, p) => c === 4 && (p === 1 || p === 2) },
  { k: 'estimator', l: 'Background estimator', t: 'select', sec: 'analysis', d: 'median', o: ['median', 'mean'], when: (c, p) => c === 4 && p === 3 },
  { k: 'bg_flank', l: 'Background flank (bp)', t: 'number', sec: 'analysis', d: 5000, min: 20, step: 500, when: (c, p) => c === 4 && p === 3 },
  { k: 'pixels', l: 'Resolution (bins)', t: 'number', sec: 'analysis', d: 800, min: 50, max: 5000, step: 50, when: C(4) },
  // Track layout
  { k: 'scale', l: 'Y scale', t: 'select', sec: 'tracks', d: 'auto_highest', o: [['auto_highest', 'Shared across tracks'], ['per_track', 'Independent per track']], when: C(0) },
  { k: 'transcripts', l: 'Transcripts', t: 'select', sec: 'tracks', d: 'collapsed', o: [['collapsed', 'Collapse per gene'], ['individual', 'Show all']], when: C(0) },
  { k: 'gene_label', l: 'Gene labels', t: 'select', sec: 'tracks', d: 'name', o: [['name', 'Gene name'], ['id', 'Gene ID'], ['none', 'None']], when: C(0) },
  { k: 'x_tick_unit', l: 'Coordinate unit', t: 'select', sec: 'tracks', d: 'auto', o: ['auto', 'bp', 'kb', 'Mb'], when: C(0) },
  { k: 'track_gap', l: 'Track spacing', t: 'range', sec: 'tracks', d: 0.18, min: 0, max: 0.8, step: 0.02, when: C(0) },
  { k: 'y_ticks', l: 'Y ticks per track', t: 'number', sec: 'tracks', d: 3, min: 2, max: 12, step: 1, when: C(0) },
  { k: 'y_decimals', l: 'Y tick decimals', t: 'number', sec: 'tracks', d: 0, min: 0, max: 6, step: 1, when: C(0) },
  // Colours & style
  { k: 'palette', l: 'Sample palette', t: 'palette', sec: 'color', d: 'studio', wide: 1 },
  { k: 'custom_palette', l: 'Custom palette', t: 'text', sec: 'color', d: '', ph: '#264653, #2a9d8f, #e9c46a', help: 'Comma-separated colours; overrides the palette above. Click a swatch in the Data panel to colour one file.', wide: 1 },
  { k: 'cmap', l: 'Heatmap colour scale', t: 'cmap', sec: 'color', d: 'auto', wide: 1, help: '“auto” fades each sample from white to its own colour. You can also type any Matplotlib name or a list of stops, e.g. #ffffff,#fdae61,#7f0000.', when: heatWhen },
  { k: 'corr_cmap', l: 'Correlation colour scale', t: 'cmap', diverging: 1, sec: 'color', d: 'RdBu_r', wide: 1, when: (c, p) => c === 2 && p === 0 },
  { k: 'pc_cmap', l: 'p-value colour scale', t: 'cmap', sec: 'color', d: 'viridis', wide: 1, when: (c, p) => c === 3 && p === 1 },
  { k: 'up_color', l: 'Up-regulated', t: 'color', sec: 'color', d: '', when: (c, p) => c === 3 && p === 0 },
  { k: 'down_color', l: 'Down-regulated', t: 'color', sec: 'color', d: '', when: (c, p) => c === 3 && p === 0 },
  { k: 'derived_color', l: 'Derived track', t: 'color', sec: 'color', d: '#3b4252', when: C(4) },
  { k: 'theme', l: 'Axis style', t: 'select', sec: 'color', d: 'clean', o: [['clean', 'Clean (open axes)'], ['classic', 'Boxed'], ['grid', 'Light grid'], ['minimal', 'Minimal']] },
  { k: 'linewidth', l: 'Line width', t: 'range', sec: 'color', d: c => (c === 0 ? 0.6 : 1.5), min: 0.2, max: 5, step: 0.1 },
  { k: 'fill_alpha', l: 'Fill opacity', t: 'range', sec: 'color', d: c => (c === 0 ? 0.85 : c === 4 ? 0.25 : 0.18), min: 0, max: 1, step: 0.01, when: c => c !== 2 && c !== 3 },
  // Text
  { k: 'title', l: 'Title', t: 'text', sec: 'text', d: () => plotName(), wide: 1 },
  { k: 'font_family', l: 'Font', t: 'combo', sec: 'text', d: () => defaultFont(), list: () => S.meta.fonts, wide: 1, help: 'Any installed font; text stays editable in SVG/PDF.' },
  { k: 'xlab', l: 'X-axis label', t: 'text', sec: 'text', d: '', ph: 'automatic', when: c => c !== 2 && c !== 3 },
  { k: 'ylab', l: 'Y-axis label', t: 'text', sec: 'text', d: '', ph: 'automatic', when: C(1) },
  { k: 'colorbar_label', l: 'Colour key label', t: 'text', sec: 'text', d: 'signal', when: heatWhen },
  { k: 'title_size', l: 'Title size', t: 'range', sec: 'text', d: 12, min: 6, max: 32, step: 0.5 },
  { k: 'axis_label_size', l: 'Axis label size', t: 'range', sec: 'text', d: 10, min: 5, max: 28, step: 0.5 },
  { k: 'tick_size', l: 'Tick label size', t: 'range', sec: 'text', d: 8, min: 4, max: 24, step: 0.5 },
  { k: 'legend_size', l: 'Legend size', t: 'range', sec: 'text', d: 9, min: 4, max: 24, step: 0.5, when: C(1) },
  { k: 'sample_label_size', l: 'Sample title size', t: 'range', sec: 'text', d: 10, min: 4, max: 28, step: 0.5, when: heatWhen },
  { k: 'track_label_size', l: 'Track label size', t: 'range', sec: 'text', d: 9, min: 4, max: 24, step: 0.5, when: C(0) },
  // Figure
  { k: 'size_unit', l: 'Units', t: 'select', sec: 'figure', d: 'mm', o: [['mm', 'millimetres'], ['in', 'inches']] },
  { k: 'exact_size', l: 'Exact size (no cropping)', t: 'check', sec: 'figure', d: false, help: 'Export at exactly width × height, e.g. 85 mm for a single journal column. Needs both values.' },
  { k: 'fig_width', l: 'Width', t: 'number', sec: 'figure', d: '', ph: 'auto', min: 1, max: 1000, step: 0.1, unitLabel: true },
  { k: 'fig_height', l: 'Height', t: 'number', sec: 'figure', d: '', ph: 'auto', min: 1, max: 1000, step: 0.1, unitLabel: true },
  { k: '_sizes', t: 'sizes', sec: 'figure', wide: 1 },
  { k: 'dpi', l: 'PNG resolution (dpi)', t: 'combo', sec: 'figure', d: 300, list: ['150', '300', '600', '1200'] },
  { k: 'transparent', l: 'Transparent background', t: 'check', sec: 'figure', d: false },
  // Advanced
  { k: '_advanced', t: 'json', sec: 'advanced', wide: 1 },
];
const FIELD = Object.fromEntries(FIELDS.map(f => [f.k, f]));
// Changing these keys can show or hide other fields, so the panel is re-rendered.
const DEP_KEYS = new Set(['rel_mode', 'rel_scope', 'alignment', 'cluster', 'autoscale', 'ln_h', 'ylim_manual', 'legend', 'sort', 'palette', 'custom_palette', 'size_unit', 'exact_size']);
// [label, unit, width, height]
const SIZE_PRESETS = [['Single column', 'mm', 85, ''], ['1.5 column', 'mm', 114, ''], ['Double column', 'mm', 180, ''], ['Slide 16:9', 'in', 11, 6.2]];
const REGION_KEYS = ['chrom', 'start', 'end', 'assembly', 'numerator', 'denominator', 'highlights'];

function defaultFont() { const f = S.meta?.fonts || []; return ['Arial', 'Helvetica', 'DejaVu Sans'].find(x => f.includes(x)) || 'DejaVu Sans'; }
function fieldDefault(f) { return typeof f.d === 'function' ? f.d(S.category, S.plot) : f.d; }
function val(k) { const f = FIELD[k]; const v = S.settings[k]; return v === undefined || v === null ? (f ? fieldDefault(f) : undefined) : v; }
function applies(f) { return !f.when || f.when(S.category, S.plot); }
function isChanged(f) { return S.settings[f.k] !== undefined && String(S.settings[f.k]) !== String(fieldDefault(f)); }
function effectiveSettings() {
  const out = {};
  for (const f of FIELDS) if (!f.k.startsWith('_')) out[f.k] = val(f.k);
  // Keep keys from loaded configs, templates or project defaults that this UI does not expose.
  for (const [k, v] of Object.entries(S.settings)) if (!(k in FIELD)) out[k] = v;
  out.flank = S.category === 4 && S.plot === 3 ? val('bg_flank') : val('flank');
  return out;
}
function setSetting(k, v, { rerender } = {}) {
  S.settings[k] = v;
  if (FIELD[k]) {
    const el = $(`.field[data-key="${k}"]`);
    if (el) el.classList.toggle('changed', isChanged(FIELD[k]));
  }
  if (k === 'palette' || k === 'custom_palette') renderFiles();
  if (rerender || DEP_KEYS.has(k)) renderSettings();
  if (k === 'title') renderCrumbs();
  if (k === 'rel_mode') { renderNeeds(); renderPlots(); }
  changed();
}
function chromList() {
  const counts = new Map();
  for (const f of enabled()) for (const c of f.info?.chroms || []) counts.set(c, (counts.get(c) || 0) + 1);
  return [...counts.entries()].sort((a, b) => b[1] - a[1]).map(([c]) => c).slice(0, 200);
}
function groupList() { return [...new Set(byRole('signal').map(f => f.group).filter(Boolean))]; }

/* ------------------------------------------------------------------ icons */
const ICONS = [
  ['<path d="M2 7c3-5 4-5 6 0s4 2 6 0 5-4 7 0 4 2 7 0" fill="none" stroke="#2a7f78" stroke-width="1.4"/><path d="M2 14c3-3 4-3 6 0s4 1 6-3 5 2 7 2 4-2 7-1" fill="none" stroke="#c0714f" stroke-width="1.4"/><rect x="2" y="17" width="26" height="1" fill="#8a95a3"/><rect x="8" y="16" width="5" height="3" fill="#4f6d8f"/>',
    '<path d="M2 16c3-9 5-9 7-3s4 4 6-3 5-8 7 0 4 4 6 2" fill="none" stroke="#2a7f78" stroke-width="1.5"/><path d="M2 17c3-5 5-6 7-1s4 2 6-5 5-3 7 1 4 2 6 1" fill="none" stroke="#c0714f" stroke-width="1.5"/>',
    '<rect x="2" y="7.5" width="26" height="1" fill="#8a95a3"/><rect x="5" y="5" width="4" height="6" fill="#4f6d8f"/><rect x="12" y="5" width="6" height="6" fill="#4f6d8f"/><rect x="21" y="5" width="3" height="6" fill="#4f6d8f"/><rect x="2" y="15.5" width="20" height="1" fill="#8a95a3"/><rect x="6" y="13" width="5" height="6" fill="#4f6d8f" opacity=".6"/><rect x="15" y="13" width="4" height="6" fill="#4f6d8f" opacity=".6"/>'],
  ['<defs><linearGradient id="g1" x1="0" x2="0" y1="0" y2="1"><stop offset="0" stop-color="#14403c"/><stop offset="1" stop-color="#e5f1ef"/></linearGradient></defs><rect x="4" y="3" width="10" height="16" fill="url(#g1)"/><rect x="16" y="3" width="10" height="16" fill="url(#g1)" opacity=".7"/>',
    '<path d="M2 18c4 0 6-1 8-6s3-9 5-9 3 4 5 9 4 6 8 6" fill="none" stroke="#2a7f78" stroke-width="1.6"/><path d="M2 18c4 0 7 0 9-3s3-5 4-5 2 2 4 5 5 3 9 3" fill="none" stroke="#c0714f" stroke-width="1.4"/>',
    '<path d="M4 7c2 0 3-4 5-4s3 4 5 4" fill="none" stroke="#2a7f78" stroke-width="1.3"/><path d="M16 7c2 0 3-3 5-3s3 3 5 3" fill="none" stroke="#c0714f" stroke-width="1.3"/><rect x="4" y="9" width="10" height="11" fill="url(#g1)"/><rect x="16" y="9" width="10" height="11" fill="#c0714f" opacity=".55"/>',
    '<rect x="4" y="3" width="20" height="5" fill="#14403c"/><rect x="4" y="9" width="20" height="5" fill="#2a7f78"/><rect x="4" y="15" width="20" height="4" fill="#9fd0c8"/><rect x="25" y="3" width="2" height="5" fill="#c0714f"/><rect x="25" y="9" width="2" height="5" fill="#4f6d8f"/><rect x="25" y="15" width="2" height="4" fill="#d9a441"/>'],
  ['<rect x="6" y="2" width="6" height="6" fill="#b2182b"/><rect x="12" y="2" width="6" height="6" fill="#f4a582"/><rect x="18" y="2" width="6" height="6" fill="#92c5de"/><rect x="6" y="8" width="6" height="6" fill="#f4a582"/><rect x="12" y="8" width="6" height="6" fill="#b2182b"/><rect x="18" y="8" width="6" height="6" fill="#d1e5f0"/><rect x="6" y="14" width="6" height="6" fill="#92c5de"/><rect x="12" y="14" width="6" height="6" fill="#d1e5f0"/><rect x="18" y="14" width="6" height="6" fill="#b2182b"/>',
    '<path d="M4 19H27M4 19V2" stroke="#8a95a3" stroke-width=".8"/><circle cx="9" cy="14" r="2" fill="#2a7f78"/><circle cx="12" cy="12" r="2" fill="#2a7f78"/><circle cx="20" cy="6" r="2" fill="#c0714f"/><circle cx="23" cy="8" r="2" fill="#c0714f"/><circle cx="16" cy="15" r="2" fill="#4f6d8f"/>',
    '<rect x="3" y="2" width="11" height="8" fill="none" stroke="#8a95a3" stroke-width=".6"/><rect x="16" y="2" width="11" height="8" fill="none" stroke="#8a95a3" stroke-width=".6"/><rect x="3" y="12" width="11" height="8" fill="none" stroke="#8a95a3" stroke-width=".6"/><rect x="16" y="12" width="11" height="8" fill="none" stroke="#8a95a3" stroke-width=".6"/><path d="M5 9l2-2 2-1 3-3" stroke="#2a7f78" stroke-width="1.2" fill="none" stroke-dasharray="1 1"/><path d="M18 18l3-3 3-1 2-1" stroke="#2a7f78" stroke-width="1.2" fill="none" stroke-dasharray="1 1"/><rect x="5" y="15" width="2" height="5" fill="#9fd0c8"/><rect x="8" y="13" width="2" height="7" fill="#9fd0c8"/><rect x="11" y="16" width="2" height="4" fill="#9fd0c8"/>'],
  ['<circle cx="15" cy="17" r="1.3" fill="#b9c0c8"/><circle cx="13" cy="15" r="1.3" fill="#b9c0c8"/><circle cx="17" cy="15.5" r="1.3" fill="#b9c0c8"/><circle cx="14.5" cy="13" r="1.3" fill="#b9c0c8"/><circle cx="8" cy="9" r="1.5" fill="#2a7f78"/><circle cx="5" cy="5" r="1.5" fill="#2a7f78"/><circle cx="22" cy="8" r="1.5" fill="#c0714f"/><circle cx="25" cy="4" r="1.5" fill="#c0714f"/><circle cx="20" cy="11" r="1.5" fill="#c0714f"/>',
    '<rect x="5" y="2" width="6" height="6" fill="#fde725"/><rect x="11" y="2" width="6" height="6" fill="#21918c"/><rect x="17" y="2" width="6" height="6" fill="#3b528b"/><rect x="11" y="8" width="6" height="6" fill="#fde725"/><rect x="17" y="8" width="6" height="6" fill="#5ec962"/><rect x="17" y="14" width="6" height="6" fill="#fde725"/>'],
  ['<path d="M2 8c3-5 5-5 7 0s4 2 6-1 5-4 7 0 4 2 6 1" fill="none" stroke="#2a7f78" stroke-width="1.3"/><path d="M2 9c3-2 5-3 7 0s4 0 6-3 5 0 7 1 4 0 6 0" fill="none" stroke="#c0714f" stroke-width="1.3"/><path d="M2 16h26" stroke="#8a95a3" stroke-width=".6" stroke-dasharray="1.5 1.5"/><path d="M4 16v-3M8 16v-5M12 16v2M16 16v-4M20 16v3M24 16v-2" stroke="#3b4252" stroke-width="2"/>',
    '<path d="M2 11h26" stroke="#8a95a3" stroke-width=".6" stroke-dasharray="1.5 1.5"/><path d="M2 11c2 0 3-6 5-6s3 10 5 10 3-8 5-8 3 6 5 6 4-4 6-3" fill="none" stroke="#3b4252" stroke-width="1.5"/>',
    '<path d="M2 11h26" stroke="#8a95a3" stroke-width=".6" stroke-dasharray="1.5 1.5"/><path d="M2 11c2 0 3-8 5-8s3 14 5 14 3-11 5-11 3 8 5 8 4-5 6-4" fill="none" stroke="#3b4252" stroke-width="1.5"/><text x="21" y="7" font-size="6" fill="#c0714f" font-family="sans-serif">log</text>',
    '<path d="M2 14h26" stroke="#c0714f" stroke-width=".8" stroke-dasharray="2 1.5"/><path d="M2 16c3-1 4-12 7-12s3 9 6 9 4-6 6-6 4 6 7 7" fill="#2a7f78" fill-opacity=".25" stroke="#2a7f78" stroke-width="1.4"/>'],
];
const icon = (c, p) => `<svg viewBox="0 0 30 22" aria-hidden="true">${ICONS[c][p]}</svg>`;

/* ------------------------------------------------------------------ rendering: header */
function renderTopbar() {
  const sel = $('#project');
  fill(sel, ...S.projects.map(p => h('option', { value: p.id, selected: p.id === S.project?.id },
    `${p.demo ? '◇ ' : ''}${p.name}${p.role === 'viewer' ? ' (view only)' : ''}`)));
  sel.title = S.project ? `${S.project.name}${S.project.description ? ' · ' + S.project.description : ''}` : '';
  const me = S.me.user;
  const initials = (me.name || me.email).split(/[\s.@_-]+/).filter(Boolean).slice(0, 2).map(s => s[0].toUpperCase()).join('');
  $('#userbtn').textContent = initials || '?';
  fill($('#usermenu'), 
    h('div', { class: 'label' }, me.email === 'local' ? 'Local single-user mode' : me.email),
    me.is_admin ? h('a', { href: '/admin' }, 'Admin console', h('span', { class: 'sub' }, 'users · projects · usage')) : null,
    S.me.auth === 'local' ? h('button', { onclick: async () => { await api('/api/logout', { method: 'POST', body: '{}' }); location.href = '/login'; } }, 'Sign out') : null,
  );
  $('#rolechip').hidden = !readOnly();
  $('#dropzone').hidden = readOnly();
  $('#live').disabled = readOnly();
}
function renderCategories() {
  fill($('#categories'), ...S.meta.catalog.map(([name], i) => h('button', {
    'aria-pressed': String(S.category === i),
    onclick: () => { if (S.category === i) return; S.category = i; S.plot = i === 1 ? 2 : 0; afterPlotChange(); },
  }, name)));
}
function renderPlots() {
  const plots = S.meta.catalog[S.category][1];
  fill($('#plots'), ...plots.map((name, p) => {
    const ok = requirement(S.category, p).ok;
    return h('button', {
      class: 'plot-card', role: 'tab', 'aria-selected': String(S.plot === p),
      title: ok ? name : `${name}: ${needText(S.category, p)}`,
      onclick: () => { if (S.plot === p) return; S.plot = p; afterPlotChange(); },
    }, h('span', { html: icon(S.category, p) }).firstChild, name, h('i', { class: 'ok' + (ok ? '' : ' no'), title: ok ? 'Inputs ready' : 'More inputs needed' }));
  }));
}
function needText(c, p) {
  return 'needs ' + requirement(c, p).items.map(i => `${i.n} ${ROLE_NAMES[i.role][i.n > 1 ? 1 : 0]}`).join(' + ');
}
function renderCrumbs() {
  $('#crumb-cat').textContent = S.meta.catalog[S.category][0];
  $('#crumb-plot').textContent = val('title') || plotName();
  const fig = S.figures.find(f => f.id === S.figureId);
  $('#savedname').hidden = !fig;
  if (fig) $('#savedname').textContent = '★ ' + fig.name;
}
function afterPlotChange() {
  // Titles that were left at a plot's default name follow the newly selected plot.
  if (S.settings.title && S.meta.catalog.some(([, ps]) => ps.includes(S.settings.title))) delete S.settings.title;
  renderCategories(); renderPlots(); renderNeeds(); renderSettings(); renderCrumbs();
  changed();
}

/* ------------------------------------------------------------------ rendering: data panel */
function renderProjectInfo() {
  const p = S.project;
  const quota = p.quota ? ` · uploads ${fmtBytes(p.used || 0)} / ${fmtBytes(p.quota)}` : '';
  fill($('#projectinfo'), 
    h('div', { class: 'pi-name' }, p.name, p.genome ? h('span', { class: 'chip' }, p.genome) : null),
    p.description ? h('div', { class: 'pi-desc' }, p.description) : null,
    h('div', { class: 'pi-meta' }, `${S.files.filter(f => !f.shared && f.source !== 'upload').length} delivered files${quota}${p.expires ? ' · access until ' + p.expires : ''}`));
}
function renderNeeds() {
  const req = requirement();
  const chips = req.items.map(i => h('span', { class: 'chip ' + (i.ok ? 'ok' : 'miss') },
    i.ok ? '✓' : '•', ` ${i.have}/${i.n} ${ROLE_NAMES[i.role][i.n > 1 ? 1 : 0]}`));
  const warnings = [];
  if (S.category === 0 || S.category === 4) {
    const chrom = normChrom(val('chrom'));
    const missing = byRole('signal').filter(f => f.info?.chroms?.length && !f.info.chroms.some(c => normChrom(c) === chrom));
    if (missing.length) warnings.push(h('span', { class: 'chip err', title: missing.map(f => f.label).join(', ') }, `${val('chrom')} absent in ${missing.length} signal${missing.length > 1 ? 's' : ''}`));
  }
  fill($('#needs'), h('span', { class: 'chip' }, plotName() + ':'), ...chips, ...warnings);
  const run = $('#run');
  run.disabled = readOnly() || (!S.busy && !req.ok);
  run.title = readOnly() ? 'Your access to this project is view-only' : req.ok ? 'Plot (⌘/Ctrl + Enter)' : `${plotName()} ${needText(S.category, S.plot)}`;
}

function renderFiles() {
  const list = $('#filelist');
  const roles = [['signal', 'Signals'], ['regions', 'Regions'], ['genes', 'Annotation']];
  const kids = [];
  if (!S.files.length) kids.push(h('p', { class: 'empty' }, readOnly() ? 'No files in this project.' : 'No files yet. Add your own with the drop zone above.'));
  for (const [role, title] of roles) {
    const files = S.files.filter(f => f.role === role);
    if (!files.length) continue;
    const on = files.filter(f => f.enabled !== false).length;
    kids.push(h('div', { class: 'group-head' },
      h('span', {}, `${title} · ${on}/${files.length}`),
      h('button', { class: 'btn ghost', onclick: () => { const all = on < files.length; files.forEach(f => { f.enabled = all; }); filesChanged(); } }, on < files.length ? 'All' : 'None')));
    for (const f of files) kids.push(fileRow(f));
  }
  fill(list, ...kids);
}
function fileRow(f) {
  const color = fileColor(f);
  const open = S.openFile === f.id;
  const swatch = h('label', { class: 'swatch' + (color ? '' : ' none'), title: f.role === 'genes' ? 'Gene annotation' : 'Colour (click to change)', style: color ? { background: color } : null },
    f.role !== 'genes' ? h('input', { type: 'color', value: color || '#888888', 'aria-label': 'Colour for ' + f.label,
      oninput: e => { f.color = e.target.value; swatch.style.background = f.color; },
      onchange: () => { filesChanged(); } }) : null);
  const label = h('input', { class: 'file-label', value: f.label, 'aria-label': 'Label', title: f.name,
    onchange: e => { f.label = e.target.value.trim() || f.base.label || f.name; filesChanged({ keepRows: true }); },
    onkeydown: e => { if (e.key === 'Enter') e.target.blur(); } });
  const badge = f.source === 'upload' ? 'upload' : f.shared ? 'shared' : null;
  const row = h('div', { class: 'file' + (f.enabled === false ? ' off' : '') + (open ? ' open' : ''), draggable: 'true', 'data-id': f.id },
    h('div', { class: 'file-main' },
      h('input', { type: 'checkbox', checked: f.enabled !== false, 'aria-label': 'Use ' + f.label, onchange: e => { f.enabled = e.target.checked; filesChanged(); } }),
      swatch, label,
      h('span', { class: 'tag', title: 'Group: ' + (f.group || '') }, badge ? h('em', {}, badge + ' ') : null, f.group || ''),
      h('button', { class: 'more', title: 'Details', 'aria-expanded': String(open), onclick: () => { S.openFile = open ? null : f.id; renderFiles(); } }, open ? '×' : '⋯')),
    open ? fileDetail(f) : null);
  // Drag to reorder within a role: order sets panel order and A/B for transformations.
  row.addEventListener('dragstart', e => { e.dataTransfer.setData('text/x-gc-file', f.id); e.dataTransfer.effectAllowed = 'move'; });
  row.addEventListener('dragover', e => { if ([...e.dataTransfer.types].includes('text/x-gc-file')) { e.preventDefault(); row.classList.add('dragover'); } });
  row.addEventListener('dragleave', () => row.classList.remove('dragover'));
  row.addEventListener('drop', e => {
    const id = e.dataTransfer.getData('text/x-gc-file');
    if (!id) return;
    e.preventDefault(); e.stopPropagation(); row.classList.remove('dragover');
    const from = S.files.findIndex(x => x.id === id);
    if (from < 0 || id === f.id) return;
    const [moved] = S.files.splice(from, 1);
    S.files.splice(S.files.findIndex(x => x.id === f.id), 0, moved);
    filesChanged();
  });
  return row;
}
function fileDetail(f) {
  const info = f.info || {};
  const meta = [info.kind, info.rows != null ? `${fmtBp(info.rows)} ${info.kind === 'GTF' ? 'transcripts' : 'rows'}` : null,
    info.size != null ? fmtBytes(info.size) : null, info.chroms ? `${info.n_chroms || info.chroms.length} chrom${(info.n_chroms || info.chroms.length) > 1 ? 's' : ''}` : null].filter(Boolean).join(' · ');
  const groups = [...new Set(S.files.map(x => x.group).filter(Boolean))];
  const gid = 'groups-' + f.id;
  const edited = EDITABLE.some(k => (f[k] || '') !== (f.base[k] || ''));
  return h('div', { class: 'file-detail' },
    h('div', { class: 'row' }, h('span', {}, 'Role'),
      h('select', { class: 'input', onchange: e => { f.role = e.target.value; S.openFile = f.id; filesChanged(); } },
        ...[['signal', 'Signal (BigWig / bedGraph)'], ['regions', 'Regions (BED)'], ['genes', 'Genes (GTF)']].map(([v, t]) => h('option', { value: v, selected: f.role === v }, t)))),
    h('div', { class: 'row' }, h('span', {}, 'Group'),
      h('input', { class: 'input', value: f.group || '', list: gid, placeholder: 'condition / sample group', onchange: e => { f.group = e.target.value.trim(); filesChanged(); } }),
      h('datalist', { id: gid }, ...groups.map(g => h('option', { value: g })))),
    f.role === 'signal' ? h('div', { class: 'row' }, h('span', {}, 'Heatmap'),
      h('input', { class: 'input', value: f.cmap || '', list: 'cmap-list', placeholder: 'use global colour scale', onchange: e => { f.cmap = e.target.value.trim(); filesChanged({ keepRows: true }); } })) : null,
    h('div', { class: 'meta' }, h('div', {}, f.name), meta,
      info.locus ? h('div', {}, `data: ${info.locus.chrom}:${fmtBp(info.locus.start)}–${fmtBp(info.locus.end)}`) : null,
      info.chroms?.length ? h('div', {}, 'chroms: ' + info.chroms.slice(0, 8).join(', ') + (info.chroms.length > 8 ? ' …' : '')) : null,
      f.source === 'manifest' ? h('div', {}, 'Delivered by the core (labels from samples.tsv)') : null),
    info.error ? h('div', { class: 'warn' }, info.error) : null,
    h('div', { class: 'actions' },
      edited ? h('button', { class: 'btn small ghost', onclick: () => { for (const k of EDITABLE) f[k] = f.base[k]; filesChanged(); } }, 'Reset') : null,
      info.locus ? h('button', { class: 'btn small', onclick: () => goTo(info.locus.chrom, info.locus.start, info.locus.end) }, 'Go to data') : null,
      f.deletable ? h('button', { class: 'btn small danger', onclick: () => removeFile(f) }, 'Delete upload') : null));
}
function filesChanged({ keepRows } = {}) {
  if (!keepRows) renderFiles();
  renderNeeds(); renderPlots();
  if (S.category === 3 || S.category === 0 || S.category === 4) renderSettings();
  changed();
}
async function removeFile(f) {
  try { await papi('/files/' + f.id, { method: 'DELETE' }); } catch (e) { toast(e.message); return; }
  S.files = S.files.filter(x => x.id !== f.id);
  S.openFile = null;
  filesChanged();
  toast(`Deleted ${f.label}`);
}
function hydrate(file, overrides = {}) {
  const base = { role: file.role, label: file.label, group: file.group, color: file.color, cmap: file.cmap };
  return { ...file, base, ...base, ...Object.fromEntries(EDITABLE.filter(k => k in overrides).map(k => [k, overrides[k]])), enabled: overrides.enabled ?? !file.shared };
}

/* ------------------------------------------------------------------ uploads */
function uploadFile(file) {
  return new Promise((resolve, reject) => {
    const row = h('div', { class: 'upload-row' }, h('span', {}, file.name), h('span', { class: 'muted' }, fmtBytes(file.size)), h('div', { class: 'bar' }, h('i')));
    $('#uploads').append(row);
    const xhr = new XMLHttpRequest();
    xhr.open('POST', `/api/projects/${encodeURIComponent(S.project.id)}/uploads`);
    xhr.setRequestHeader('X-GenomeCanvas', 'studio');
    xhr.setRequestHeader('X-Filename', encodeURIComponent(file.name));
    xhr.upload.onprogress = e => { if (e.lengthComputable) row.querySelector('i').style.width = (100 * e.loaded / e.total).toFixed(1) + '%'; };
    xhr.onload = () => {
      let data = {};
      try { data = JSON.parse(xhr.responseText); } catch { /* non-JSON error */ }
      if (xhr.status >= 200 && xhr.status < 300) { row.remove(); resolve(data); return; }
      row.classList.add('error'); row.lastChild.remove(); row.append(h('span', { style: { gridColumn: '1/-1' } }, data.error || 'Upload failed'));
      setTimeout(() => row.remove(), 9000);
      reject(new Error(data.error || 'Upload failed'));
    };
    xhr.onerror = () => { row.remove(); reject(new Error('Network error during upload')); };
    xhr.send(file);
  });
}
async function uploadAll(list) {
  if (readOnly()) { toast('Your access to this project is view-only'); return; }
  const files = [...list];
  if (!files.length) return;
  const max = S.meta.max_upload;
  const added = [];
  for (const file of files) {
    if (file.size > max) { toast(`${file.name} is larger than ${fmtBytes(max)}; ask the core to add it to the project instead.`); continue; }
    try { added.push(await uploadFile(file)); } catch (e) { status(e.message, 'error'); }
  }
  if (!added.length) return;
  for (const rec of added) S.files.push(hydrate(rec, { enabled: true }));
  S.project.used = (S.project.used || 0) + added.reduce((n, r) => n + (r.info?.size || 0), 0);
  renderProjectInfo();
  filesChanged();
  status(`Added ${added.length} file${added.length > 1 ? 's' : ''}.`);
  suggestLocus(added);
}
/** Point the locus at newly added data if the current window does not overlap it. */
function suggestLocus(records) {
  if (S.category !== 0 && S.category !== 4) return;
  const loc = records.map(r => r.info?.locus).find(Boolean);
  if (loc && !chromList().some(c => normChrom(c) === normChrom(val('chrom')))) goTo(loc.chrom, loc.start, loc.end);
}

/* ------------------------------------------------------------------ saved figures */
async function loadFigures() {
  try { S.figures = (await papi('/figures')).figures; } catch (e) { S.figures = []; status(e.message, 'error'); }
  renderFigures(); renderCrumbs();
}
function renderFigures() {
  $('#figcount').textContent = S.figures.length ? S.figures.length : '';
  const box = $('#figures');
  if (!S.figures.length) {
    fill(box, h('p', { class: 'empty' }, readOnly() ? 'No figures have been saved in this project yet.' : 'Plot something, then press Save to keep it here with its settings.'));
    return;
  }
  fill(box, ...S.figures.map(fig => {
    const plot = fig.config?.category != null ? S.meta.catalog[fig.config.category]?.[1][fig.config.plot] : '';
    return h('div', { class: 'figure-card' + (fig.id === S.figureId ? ' current' : '') },
      h('button', { class: 'figure-thumb', title: 'Open with its settings', onclick: () => openFigure(fig) },
        fig.thumb ? h('img', { src: fig.thumb, alt: '', loading: 'lazy' }) : h('span', {}, 'result expired: open and plot again')),
      h('div', { class: 'figure-meta' },
        h('strong', {}, fig.name),
        h('small', {}, [plot, fig.author, fmtDate(fig.updated)].filter(Boolean).join(' · ')),
        fig.notes ? h('p', {}, fig.notes) : null),
      readOnly() ? null : h('div', { class: 'figure-actions' },
        h('button', { class: 'btn ghost small', onclick: () => editFigure(fig) }, 'Rename'),
        h('button', { class: 'btn ghost small danger', onclick: () => deleteFigure(fig) }, 'Delete')));
  }));
}
async function openFigure(fig) {
  applyConfig(fig.config, { quiet: true });
  S.figureId = fig.id;
  renderCrumbs(); renderFigures();
  if (fig.job_id) {
    try {
      const data = await api('/api/jobs/' + fig.job_id);
      if (data.status === 'complete') { showResult(data, payload(), S.revision); status(`Opened “${fig.name}”.`); return; }
    } catch { /* fall through */ }
  }
  status(`Opened “${fig.name}”. Its image has expired; press Plot to redraw it.`);
}
async function saveFigure() {
  if (!S.current) return;
  const existing = S.figures.find(f => f.id === S.figureId);
  const values = await ask({
    title: 'Save figure to project',
    fields: [
      { name: 'name', label: 'Name', value: existing?.name || val('title') || plotName(), required: true },
      { name: 'notes', label: 'Notes (optional)', type: 'textarea', value: existing?.notes || '', help: 'e.g. which panel of the manuscript this is' },
      existing ? { name: 'overwrite', label: `Replace “${existing.name}”`, type: 'check', value: true } : null,
    ].filter(Boolean),
  });
  if (!values) return;
  const body = { name: values.name, notes: values.notes, config: S.current.config, job_id: S.current.id };
  try {
    if (existing && values.overwrite) {
      await api('/api/figures/' + existing.id, { method: 'PATCH', body: JSON.stringify(body) });
    } else {
      const { figure } = await papi('/figures', { method: 'POST', body: JSON.stringify(body) });
      S.figureId = figure.id;
    }
    await loadFigures();
    toast(`Saved “${values.name}”. Everyone on the project can open it.`);
  } catch (e) { toast(e.message); }
}
async function editFigure(fig) {
  const values = await ask({ title: 'Edit figure', fields: [{ name: 'name', label: 'Name', value: fig.name, required: true }, { name: 'notes', label: 'Notes', type: 'textarea', value: fig.notes }] });
  if (!values) return;
  try { await api('/api/figures/' + fig.id, { method: 'PATCH', body: JSON.stringify(values) }); await loadFigures(); } catch (e) { toast(e.message); }
}
async function deleteFigure(fig) {
  const values = await ask({ title: `Delete “${fig.name}”?`, ok: 'Delete', html: h('p', {}, 'The saved settings and image are removed for everyone on the project.') });
  if (!values) return;
  try {
    await api('/api/figures/' + fig.id, { method: 'DELETE' });
    if (S.figureId === fig.id) S.figureId = null;
    await loadFigures();
  } catch (e) { toast(e.message); }
}

/* ------------------------------------------------------------------ templates */
function renderTemplates() {
  const menu = $('#tplmenu');
  const items = S.templates.map(t => h('button', { onclick: () => applyTemplate(t), title: t.description || t.name },
    h('span', {}, t.name), h('span', { class: 'sub' }, `${t.scope === 'project' ? 'project' : 'core'} · ${S.meta.catalog[t.config.category]?.[1][t.config.plot] || ''}`)));
  fill(menu, 
    h('div', { class: 'label' }, 'Apply a template'),
    ...(items.length ? items : [h('div', { class: 'label' }, 'None yet')]),
    S.me.user.is_admin ? h('hr') : null,
    S.me.user.is_admin ? h('button', { onclick: saveTemplate }, 'Save current settings as template…') : null);
}
function applyTemplate(t) {
  closeMenus();
  const keep = Object.fromEntries(REGION_KEYS.filter(k => k in S.settings).map(k => [k, S.settings[k]]));
  S.category = t.config.category ?? S.category;
  S.plot = t.config.plot ?? S.plot;
  S.settings = { ...keep, ...(t.config.settings || {}) };
  S.advanced = typeof t.config.advanced === 'string' ? t.config.advanced : t.config.advanced && Object.keys(t.config.advanced).length ? JSON.stringify(t.config.advanced, null, 2) : '';
  S.figureId = null;
  renderAll(); changed();
  toast(`Applied template “${t.name}”. Your data selection and region were kept.`);
}
async function saveTemplate() {
  closeMenus();
  const values = await ask({
    title: 'Save as template', fields: [
      { name: 'name', label: 'Name', value: plotName(), required: true },
      { name: 'description', label: 'Description', type: 'textarea' },
      { name: 'scope', label: 'Available in', type: 'select', value: 'global', options: [['global', 'All projects'], ['project', `This project only (${S.project.name})`]] },
    ],
  });
  if (!values) return;
  const settings = { ...S.settings };
  REGION_KEYS.forEach(k => delete settings[k]);
  try {
    await api('/api/templates', { method: 'POST', body: JSON.stringify({ name: values.name, description: values.description, project_id: values.scope === 'project' ? S.project.id : null, config: { category: S.category, plot: S.plot, settings, advanced: S.advanced } }) });
    S.templates = (await papi('')).templates;
    renderTemplates();
    toast('Template saved');
  } catch (e) { toast(e.message); }
}

/* ------------------------------------------------------------------ rendering: settings */
function matchesSearch(f, q) {
  if (!q) return true;
  if (f.k === '_advanced') return 'advanced json override rcparams matplotlib'.includes(q) || q.includes('json');
  if (f.k.startsWith('_')) return false;
  return `${f.l || ''} ${f.k} ${f.help || ''} ${f.sec}`.toLowerCase().includes(q);
}
function renderSettings() {
  const body = $('#settingsbody');
  const focusKey = document.activeElement?.closest?.('.field')?.dataset.key;
  const q = S.search.trim().toLowerCase();
  const kids = [];
  for (const [id, title] of SECTIONS) {
    const fields = FIELDS.filter(f => f.sec === id && applies(f) && matchesSearch(f, q));
    if (!fields.length) continue;
    const changedCount = fields.filter(isChanged).length;
    const det = h('details', { class: 'section', open: q || S.sections[id] ? true : null, 'data-sec': id },
      h('summary', {}, title, changedCount ? h('span', { class: 'count' }, `${changedCount} changed`) : null),
      h('div', { class: 'section-body' }, ...fields.map(renderField)));
    det.addEventListener('toggle', () => { if (!q) { S.sections[id] = det.open; persist(); } });
    kids.push(det);
  }
  if (!kids.length) kids.push(h('p', { class: 'empty' }, 'No settings match “' + S.search + '”.'));
  fill(body, ...kids, h('datalist', { id: 'cmap-list' }, ...Object.keys(S.meta.colormaps).map(n => h('option', { value: n }))));
  if (focusKey) $(`.field[data-key="${focusKey}"] :is(input,select,textarea)`)?.focus();
}

function renderField(f) {
  const v = val(f.k);
  const wrap = h('div', { class: `field${f.wide ? ' wide' : ''}${isChanged(f) ? ' changed' : ''}${f.dim && f.dim() ? ' dim' : ''}`, 'data-key': f.k });
  const id = 'f-' + f.k;
  const text = f.unitLabel ? `${f.l} (${unit()})` : f.l;
  const label = () => h('label', { for: id }, text);
  const help = f.help ? h('div', { class: 'help' }, f.help) : null;
  switch (f.t) {
    case 'check':
      wrap.classList.add('check');
      put(wrap, h('input', { type: 'checkbox', id, checked: !!(v === true || v === 'true'), onchange: e => setSetting(f.k, e.target.checked) }), label());
      if (help) put(wrap, help);
      if (f.k === 'exact_size' && v && !(val('fig_width') && val('fig_height'))) put(wrap, h('div', { class: 'field-note' }, 'Set both width and height for an exact size.'));
      return wrap;
    case 'select': {
      const sel = h('select', { id, onchange: e => setSetting(f.k, e.target.value) },
        ...(typeof f.o === 'function' ? f.o() : f.o).map(o => { const [value, t] = Array.isArray(o) ? o : [o, o]; return h('option', { value, selected: String(v) === String(value) }, t); }));
      return put(wrap, label(), sel, help);
    }
    case 'combo': {
      const list = typeof f.list === 'function' ? f.list() : f.list;
      const lid = 'dl-' + f.k;
      return put(wrap, label(), h('input', { id, value: v ?? '', list: lid, placeholder: f.ph || '', autocomplete: 'off', onchange: e => setSetting(f.k, e.target.value.trim()) }),
        h('datalist', { id: lid }, ...(list || []).map(x => h('option', { value: x }))), help);
    }
    case 'number':
    case 'text':
      return put(wrap, label(), h('input', {
        id, type: f.t, value: v ?? '', placeholder: f.ph || '', min: f.min, max: f.max, step: f.step,
        onchange: e => setSetting(f.k, f.t === 'number' && e.target.value !== '' ? Number(e.target.value) : e.target.value),
        onkeydown: e => { if (e.key === 'Enter') e.target.dispatchEvent(new Event('change')); },
      }), help);
    case 'textarea':
      return put(wrap, label(), h('textarea', { id, rows: 3, placeholder: f.ph || '', value: v ?? '', onchange: e => setSetting(f.k, e.target.value) }), help);
    case 'range': {
      const out = h('output', {}, String(v));
      return put(wrap, h('label', { for: id }, text, out), h('input', {
        id, type: 'range', min: f.min, max: f.max, step: f.step, value: v,
        oninput: e => { out.textContent = e.target.value; },
        onchange: e => setSetting(f.k, Number(e.target.value)),
      }), help);
    }
    case 'color': {
      const fallback = f.k === 'up_color' ? paletteColors()[1 % paletteColors().length] : f.k === 'down_color' ? paletteColors()[0] : '#3b4252';
      const input = h('input', { value: v || '', placeholder: v ? '' : 'palette', onchange: e => setSetting(f.k, e.target.value.trim()) });
      return put(wrap, label(), h('div', { class: 'color-input' },
        h('input', { type: 'color', id, value: /^#[0-9a-f]{6}$/i.test(v) ? v : fallback, oninput: e => { input.value = e.target.value; }, onchange: e => setSetting(f.k, e.target.value) }), input));
    }
    case 'palette': {
      const pals = Object.entries(S.meta.palettes);
      return put(wrap, h('label', {}, f.l), h('div', { class: 'palettes' }, ...pals.map(([key, p]) => h('button', {
        class: 'palette', type: 'button', 'aria-pressed': String(v === key && !val('custom_palette')),
        onclick: () => { S.settings.custom_palette = ''; setSetting('palette', key, { rerender: true }); },
      }, h('span', { class: 'name' }, p.name), h('span', { class: 'dots' }, ...p.colors.slice(0, 8).map(c => h('i', { style: { background: c } })))))));
    }
    case 'cmap': {
      const maps = f.diverging ? S.meta.diverging : S.meta.colormaps;
      const input = h('input', { id, value: v, list: 'cmap-list', placeholder: 'name or #hex,#hex,…', onchange: e => setSetting(f.k, e.target.value.trim(), { rerender: true }) });
      return put(wrap, label(), h('div', { class: 'cmaps' }, ...Object.entries(maps).map(([name, stops]) => h('button', {
        class: 'cmap', type: 'button', 'aria-pressed': String(v === name), title: name,
        onclick: () => setSetting(f.k, name, { rerender: true }),
      }, h('i', { style: { background: name === 'auto' ? autoGradient() : `linear-gradient(90deg,${stops.join(',')})` } }), h('span', {}, name)))), input, help);
    }
    case 'locus': return locusWidget(wrap);
    case 'groups': {
      const counts = {};
      byRole('signal').forEach(s => { counts[s.group || '(none)'] = (counts[s.group || '(none)'] || 0) + 1; });
      const ok = (counts[val('numerator')] || 0) >= 2 && (counts[val('denominator')] || 0) >= 2 && val('numerator') !== val('denominator');
      return put(wrap, h('div', { class: ok ? 'help' : 'field-note' },
        'Signal groups: ' + (Object.entries(counts).map(([g, n]) => `${g} ×${n}`).join(', ') || 'none') +
        (ok ? '' : '. Each side needs ≥2 replicates; set groups in the Data panel (⋯).')));
    }
    case 'sizes':
      return put(wrap, h('div', { class: 'inline-btns' }, ...SIZE_PRESETS.map(([name, u, w, hgt]) => h('button', {
        class: 'btn small', type: 'button', title: `${w}${hgt ? ' × ' + hgt : ''} ${u}`,
        onclick: () => { S.settings.size_unit = u; S.settings.fig_width = w; S.settings.fig_height = hgt; setSetting('fig_width', w, { rerender: true }); },
      }, name)), h('button', { class: 'btn small ghost', type: 'button', onclick: () => { S.settings.fig_width = ''; S.settings.exact_size = false; setSetting('fig_height', '', { rerender: true }); } }, 'Auto')));
    case 'json': return advancedWidget(wrap);
  }
  return wrap;
}
function autoGradient() {
  const c = paletteColors().slice(0, 3);
  return `linear-gradient(90deg,#fff,${c[0]} 33%,#fff 34%,${c[1] || c[0]} 66%,#fff 67%,${c[2] || c[0]})`;
}

function locusWidget(wrap) {
  const results = h('div', { class: 'gene-results', hidden: true });
  const genesFiles = byRole('genes');
  const search = h('input', {
    class: 'input', type: 'search', placeholder: genesFiles.length ? 'Gene name or chr:start-end' : 'chr:start-end (enable a GTF to search genes)',
    oninput: debounce(async e => {
      const q = e.target.value.trim();
      const loc = parseLocus(q);
      const items = [];
      if (loc) items.push(h('button', { onclick: () => { goTo(loc.chrom, loc.start, loc.end, 0); results.hidden = true; } }, `Go to ${loc.chrom}:${fmtBp(loc.start)}-${fmtBp(loc.end)}`, h('small', {}, 'locus')));
      else if (q.length >= 2 && genesFiles.length) {
        try {
          const { hits } = await papi(`/genes?files=${genesFiles.map(f => f.id).join(',')}&q=${encodeURIComponent(q)}`);
          for (const g of hits) items.push(h('button', { onclick: () => { goTo(g.chrom, g.start, g.end, 0.15); results.hidden = true; search.value = g.name; } },
            h('span', {}, g.name, g.strand ? ` (${g.strand})` : ''), h('small', {}, `${g.chrom}:${fmtBp(g.start)}`)));
          if (!hits.length) items.push(h('button', { disabled: true }, 'No matching genes'));
        } catch (err) { items.push(h('button', { disabled: true }, err.message)); }
      }
      fill(results, ...items);
      results.hidden = !items.length;
    }, 220),
    onkeydown: e => { if (e.key === 'Enter') results.querySelector('button:not([disabled])')?.click(); },
  });
  const start = Number(val('start')), end = Number(val('end')), span = Math.max(1, end - start);
  const nav = (newStart, newEnd) => goTo(val('chrom'), Math.max(0, Math.round(newStart)), Math.round(newEnd), 0);
  return put(wrap,
    h('label', {}, 'Navigate', h('output', {}, `${fmtBp(span)} bp`)), search, results,
    h('div', { class: 'inline-btns' },
      h('button', { class: 'btn small', type: 'button', title: 'Pan left', onclick: () => nav(start - span / 2, end - span / 2) }, '←'),
      h('button', { class: 'btn small', type: 'button', title: 'Zoom out 2×', onclick: () => nav(start - span / 2, end + span / 2) }, '−'),
      h('button', { class: 'btn small', type: 'button', title: 'Zoom in 2×', onclick: () => nav(start + span / 4, end - span / 4) }, '+'),
      h('button', { class: 'btn small', type: 'button', title: 'Pan right', onclick: () => nav(start + span / 2, end + span / 2) }, '→'),
      h('button', { class: 'btn small', type: 'button', title: 'Show the full extent of the enabled signals', onclick: fitToData }, 'Fit to data')));
}
function parseLocus(text) {
  const m = String(text).replace(/,/g, '').match(/^([\w.]+):(\d+)\s*[-–]\s*(\d+)$/);
  return m && +m[3] > +m[2] ? { chrom: m[1], start: +m[2], end: +m[3] } : null;
}
function goTo(chrom, start, end, pad = 0) {
  const p = Math.round(Math.max((end - start) * pad, pad ? 500 : 0));
  S.settings.chrom = chrom; S.settings.start = Math.max(0, start - p); S.settings.end = end + p;
  if (S.category !== 0 && S.category !== 4) { S.category = 0; S.plot = 0; renderCategories(); renderPlots(); }
  renderSettings(); renderNeeds(); changed();
  status(`Region set to ${chrom}:${fmtBp(S.settings.start)}–${fmtBp(S.settings.end)}`);
}
function fitToData() {
  const loci = (byRole('signal').length ? byRole('signal') : enabled()).map(f => f.info?.locus).filter(Boolean);
  if (!loci.length) { toast('No coordinate summary available for the enabled files.'); return; }
  const counts = {};
  loci.forEach(l => { counts[l.chrom] = (counts[l.chrom] || 0) + 1; });
  const chrom = Object.entries(counts).sort((a, b) => b[1] - a[1])[0][0];
  const on = loci.filter(l => l.chrom === chrom);
  goTo(chrom, Math.min(...on.map(l => l.start)), Math.max(...on.map(l => l.end)), 0.02);
}

function advancedWidget(wrap) {
  const fields = S.meta.style_fields;
  const err = h('div', { class: 'help', style: { color: 'var(--danger)' } });
  const area = h('textarea', {
    rows: 8, spellcheck: 'false', value: S.advanced,
    placeholder: '{\n  "rc": {"axes.linewidth": 0.6},\n  "heatmap": {"profile_ratio": 0.35},\n  "profile": {"figsize": [6, 3.5]},\n  "layout": {"show_cytoband": false},\n  "track_style": {"fill": true}\n}',
    oninput: e => { try { if (e.target.value.trim()) JSON.parse(e.target.value); err.textContent = ''; } catch (x) { err.textContent = x.message; } },
    onchange: e => { S.advanced = e.target.value; wrap.classList.toggle('changed', !!S.advanced.trim()); changed(); },
  });
  return put(wrap, h('label', {}, 'Overrides (JSON)'), area, err, h('div', { class: 'help' },
    'Set any engine option the panels do not expose. Sections: ', h('code', {}, 'rc'), ' (Matplotlib rcParams: axes., lines., font., …), ',
    h('code', {}, 'heatmap'), ', ', h('code', {}, 'profile'), ', ', h('code', {}, 'layout'), ' / ', h('code', {}, 'style'), ' (track figures), ',
    h('code', {}, 'track_style'), '.'),
  h('details', { class: 'help' }, h('summary', {}, 'Available heatmap / profile fields'),
    h('div', {}, h('strong', {}, 'heatmap: '), fields.heatmap.join(', ')), h('div', {}, h('strong', {}, 'profile: '), fields.profile.join(', '))));
}

/* ------------------------------------------------------------------ run & preview */
function changed() {
  S.revision++;
  persist();
  const paper = $('#paper');
  if (S.current) paper.classList.toggle('stale', S.resultRevision !== S.revision);
  if (S.live && !readOnly() && requirement().ok) scheduleRun();
  else if (S.current && S.resultRevision !== S.revision) status('Settings changed. Press Plot (⌘↵) to update the figure.');
}
const scheduleRun = debounce(() => run(), 650);

function payload() {
  let advanced = {};
  if (S.advanced.trim()) advanced = JSON.parse(S.advanced);
  const settings = effectiveSettings();
  return {
    category: S.category, plot: S.plot,
    files: enabled().map(f => ({ id: f.id, role: f.role, group: f.group, label: f.label, color: f.color || '', cmap: f.cmap || '' })),
    settings, title: settings.title || plotName(), palette: settings.palette, advanced,
  };
}
async function run() {
  if (readOnly()) return;
  if (S.busy) { S.queued = true; return; }
  const req = requirement();
  if (!req.ok) { status(`${plotName()} ${needText(S.category, S.plot)}.`, 'error'); return; }
  let request;
  try { request = payload(); } catch (e) { status('Advanced JSON is invalid: ' + e.message, 'error'); return; }
  const revision = S.revision;
  const config = configObject();
  S.busy = true;
  $('#busy').hidden = false; $('#busytext').textContent = 'Plotting ' + plotName() + '…';
  $('.run-label').textContent = 'Plotting…';
  status('Submitting job…');
  try {
    const job = await papi('/jobs', { method: 'POST', body: JSON.stringify(request) });
    let data;
    for (let i = 0; ; i++) {
      data = await api('/api/jobs/' + job.id);
      if (data.status === 'error') throw new Error(data.error);
      if (data.status === 'complete') break;
      $('#busytext').textContent = data.status === 'queued' ? `Waiting in queue${data.position ? ' (#' + data.position + ')' : ''}…` : `Plotting ${plotName()}…`;
      await new Promise(r => setTimeout(r, i < 10 ? 300 : 700));
    }
    data.config = config;
    showResult(data, request, revision);
    addHistory(data, request);
    const counts = data.counts ? ` · ${data.counts.signals} signal${data.counts.signals === 1 ? '' : 's'}` : '';
    status(`Done in ${data.elapsed}s${counts}${data.warnings?.length ? ' · ' + data.warnings.join(' ') : ''}`);
  } catch (e) {
    showError(e.message);
    status('Plot failed: ' + e.message, 'error');
  } finally {
    S.busy = false;
    $('#busy').hidden = true; $('.run-label').textContent = 'Plot';
    renderNeeds();
    if (S.queued) { S.queued = false; if (S.revision !== S.resultRevision && S.live) run(); }
  }
}
function showResult(data, request, revision) {
  S.current = { ...data, request, config: data.config || configObject() };
  S.resultRevision = revision;
  const img = $('#figure');
  img.onload = () => applyZoom();
  img.src = data.preview_png + '?t=' + Date.now();
  img.alt = request.title;
  img.hidden = false;
  const ph = $('#placeholder');
  ph.hidden = true; ph.classList.remove('error');
  $('#paper').classList.toggle('stale', revision !== S.revision);
  $('#savebtn').disabled = readOnly();
  renderExportMenu();
  const pca = data.files?.find(f => f.name === 'pca.json');
  S.current.interactive = pca ? pca.url : null;
  $('#viewtoggle').hidden = !pca;
  setView(pca ? S.view : 'static');
}

/* ------------------------------------------------------------------ interactive (Plotly) */
let plotlyLoading;
function loadPlotly() {
  if (window.Plotly) return Promise.resolve(window.Plotly);
  plotlyLoading ||= new Promise((resolve, reject) => {
    const tag = h('script', { src: '/vendor/plotly.min.js' });
    tag.onload = () => resolve(window.Plotly);
    tag.onerror = () => { plotlyLoading = null; reject(new Error('Interactive view unavailable (install plotly on the server)')); };
    document.head.append(tag);
  });
  return plotlyLoading;
}
async function setView(view) {
  const interactive = view === 'interactive' && S.current?.interactive;
  $$('#viewtoggle button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.view === (interactive ? 'interactive' : 'static'))));
  $('#figure').hidden = !!interactive || !S.current;
  $('#interactive').hidden = !interactive;
  $('#paper').classList.toggle('interactive-on', !!interactive);
  $('#zoom').style.visibility = interactive ? 'hidden' : '';
  if (!interactive) { applyZoom(); return; }
  try {
    const [Plotly, data] = await Promise.all([loadPlotly(), fetch(S.current.interactive).then(r => r.json())]);
    renderPca(Plotly, data);
  } catch (e) {
    toast(e.message);
    S.view = 'static'; setView('static');
  }
}
function renderPca(Plotly, data) {
  const v = data.variance;
  const font = { family: `${data.font || 'Arial'}, Helvetica, sans-serif`, size: 12, color: '#1f2933' };
  const traces = data.points.map(p => ({
    type: 'scatter', mode: 'markers', x: [p.pcs[0]], y: [p.pcs[1]], name: p.name, legendgroup: p.group,
    marker: { size: Math.max(8, Math.sqrt(data.point_size) * 1.6), color: p.color, symbol: p.symbol, line: { color: '#fff', width: 1 } },
    customdata: [[p.group || '—', p.pcs[2], p.file]],
    hovertemplate: `<b>${escapeHtml(p.name)}</b><br>group: %{customdata[0]}<br>PC1: %{x:.2f} (${v[0].toFixed(1)}%)<br>PC2: %{y:.2f} (${v[1].toFixed(1)}%)<br>PC3: %{customdata[1]:.2f} (${v[2].toFixed(1)}%)<br><span style="color:#8a95a3">%{customdata[2]}</span><extra></extra>`,
  }));
  const annotations = data.points.map(p => ({
    x: p.pcs[0], y: p.pcs[1], text: escapeHtml(p.name), showarrow: true, arrowwidth: 0.8, arrowcolor: '#9aa5b1',
    ax: p.label_offset[0] * 1.33 + 10, ay: -(p.label_offset[1] * 1.33 + 6), font: { ...font, size: 11 },
  }));
  const layout = {
    title: { text: escapeHtml(data.title), font: { ...font, size: 15 } }, font, showlegend: false, hovermode: 'closest',
    margin: { l: 64, r: 24, t: 56, b: 64 }, paper_bgcolor: '#fff', plot_bgcolor: '#fff', annotations, autosize: true,
    xaxis: { title: { text: `PC1 (${v[0].toFixed(1)}%)` }, zeroline: true, zerolinecolor: '#e3e7ec', showline: true, linecolor: '#52606d', ticks: 'outside', gridcolor: '#f3f4f6' },
    yaxis: { title: { text: `PC2 (${v[1].toFixed(1)}%)` }, zeroline: true, zerolinecolor: '#e3e7ec', showline: true, linecolor: '#52606d', ticks: 'outside', gridcolor: '#f3f4f6' },
  };
  const config = {
    responsive: true, displaylogo: false, editable: false, edits: { annotationTail: true }, // drag a label; its arrowhead stays on the point
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: { format: 'svg', filename: slug(data.title) },
  };
  Plotly.react($('#interactive'), traces, layout, config);
}
const escapeHtml = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
function showError(message) {
  const ph = $('#placeholder');
  if (S.current) { toast(message); return; }
  ph.hidden = false; ph.classList.add('error');
  fill(ph, h('span', { html: '<svg viewBox="0 0 120 64"><path d="M4 52c10-30 18-30 24-12s12 18 18-6 12-30 20-8 12 22 18 6 10-18 18-12 10 10 14 6"/></svg>' }).firstChild,
    h('h3', {}, 'The figure could not be drawn'), h('pre', {}, message));
}
function addHistory(data, request) {
  S.history.unshift({ data, request, snapshot: { ...configObject(), figureId: S.figureId } });
  S.history = S.history.slice(0, 16);
  renderHistory();
}
function renderHistory() {
  fill($('#history'), ...S.history.map((item, i) => h('button', {
    class: 'thumb' + (S.current?.preview_png === item.data.preview_png ? ' current' : ''), title: `${item.request.title}: click to restore these settings`,
    onclick: () => restoreSnapshot(item),
  }, h('img', { src: item.data.preview_png, alt: '', loading: 'lazy' }), h('span', {}, `${S.history.length - i}. ${S.meta.catalog[item.snapshot.category][1][item.snapshot.plot]}`))));
}
function restoreSnapshot(item) {
  applyConfig(item.snapshot, { quiet: true });
  S.figureId = item.snapshot.figureId;
  showResult(item.data, item.request, S.revision);
  renderHistory(); renderCrumbs();
  status('Restored figure and its settings from history.');
}

function renderExportMenu() {
  const menu = $('#exportmenu');
  const r = S.current;
  if (!r) { fill(menu, h('div', { class: 'label' }, 'Plot first to download')); return; }
  const base = slug(r.request.title);
  const dl = (f, name) => `${f.url}?download=${encodeURIComponent(name)}`;
  const fig = ['png', 'svg', 'pdf'].map(ext => r.files.find(f => f.name === 'figure.' + ext)).filter(Boolean);
  const repro = ['methods.txt', 'request.json', 'provenance.json'];
  const data = r.files.filter(f => !f.name.startsWith('figure.') && !/\.(png|pdf|svg|yaml|html)$/.test(f.name) && !repro.includes(f.name) && f.name !== 'pca.json');
  const extra = r.files.filter(f => /\.(png|pdf|svg)$/.test(f.name) && !f.name.startsWith('figure.'));
  const link = (f, name, sub) => h('a', { href: dl(f, name), download: name }, name, h('span', { class: 'sub' }, sub ?? fmtBytes(f.size || 0)));
  const methods = r.files.find(f => f.name === 'methods.txt');
  fill(menu, 
    h('div', { class: 'label' }, 'Figure'),
    ...fig.map(f => link(f, `${base}.${f.name.split('.').pop()}`, f.name.endsWith('png') ? `${r.request.settings.dpi} dpi · ${fmtBytes(f.size || 0)}` : 'vector, editable text')),
    ...r.files.filter(f => f.name.endsWith('.html')).map(f => link(f, `${base}.interactive.html`, 'interactive · hover & drag')),
    navigator.clipboard && window.ClipboardItem ? h('button', { onclick: copyPng }, 'Copy PNG to clipboard') : null,
    extra.length ? h('div', { class: 'label' }, 'Additional panels') : null, ...extra.map(f => link(f, f.name)),
    data.length ? h('div', { class: 'label' }, 'Data') : null, ...data.map(f => link(f, f.name)),
    h('hr'), h('div', { class: 'label' }, 'Methods & reproducibility'),
    methods ? h('button', { onclick: () => copyText(methods.url, 'Methods text copied') }, 'Copy methods paragraph', h('span', { class: 'sub' }, 'for your manuscript')) : null,
    ...r.files.filter(f => repro.includes(f.name)).map(f => link(f, f.name)));
}
async function copyPng() {
  try {
    const blob = await (await fetch(S.current.preview_png)).blob();
    await navigator.clipboard.write([new ClipboardItem({ 'image/png': blob })]);
    toast('PNG copied to clipboard');
  } catch (e) { toast('Copy failed: ' + e.message); }
}
async function copyText(url, message) {
  try { await navigator.clipboard.writeText(await (await fetch(url)).text()); toast(message); } catch (e) { toast('Copy failed: ' + e.message); }
}

/* ------------------------------------------------------------------ zoom */
function applyZoom() {
  const img = $('#figure'), paper = $('#paper');
  $$('#zoom button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.zoom === 'fit' ? S.zoom === 'fit' : b.dataset.zoom === '1' && S.zoom === 1)));
  paper.classList.toggle('fit', S.zoom === 'fit');
  if (S.zoom === 'fit' || !img.naturalWidth) { img.style.width = ''; return; }
  // 1:1 means printed size on a 96-dpi screen.
  const dpi = Number(S.current?.request?.settings?.dpi) || 300;
  img.style.width = Math.round(img.naturalWidth * (96 / dpi) * S.zoom) + 'px';
}
function setZoom(z) {
  const current = S.zoom === 'fit' ? fitScale() : S.zoom;
  S.zoom = z === 'fit' ? 'fit' : z === 'in' ? Math.min(8, current * 1.25) : z === 'out' ? Math.max(0.1, current / 1.25) : Number(z);
  applyZoom();
}
function fitScale() {
  const img = $('#figure');
  const dpi = Number(S.current?.request?.settings?.dpi) || 300;
  return img.naturalWidth ? img.clientWidth / (img.naturalWidth * 96 / dpi) : 1;
}

/* ------------------------------------------------------------------ config import/export */
function configObject() {
  return {
    app: 'GenomeCanvas Studio', version: 3, project: S.project?.id, category: S.category, plot: S.plot, plot_name: plotName(),
    settings: { ...S.settings }, advanced: S.advanced, files: S.files.map(f => ({ id: f.id, name: f.name, enabled: f.enabled !== false, ...Object.fromEntries(EDITABLE.map(k => [k, f[k] || ''])) })),
  };
}
function applyConfig(cfg, { quiet } = {}) {
  if (cfg.category != null) S.category = Number(cfg.category);
  if (cfg.plot != null) S.plot = Number(cfg.plot);
  S.settings = { ...(cfg.settings || {}) };
  if (cfg.palette) S.settings.palette = cfg.palette;
  if (cfg.title) S.settings.title = cfg.title;
  S.advanced = typeof cfg.advanced === 'string' ? cfg.advanced : cfg.advanced && Object.keys(cfg.advanced).length ? JSON.stringify(cfg.advanced, null, 2) : '';
  let missing = 0;
  if (Array.isArray(cfg.files)) {
    const order = [];
    for (const saved of cfg.files) {
      const f = S.files.find(x => x.id === saved.id) || S.files.find(x => x.name === saved.name);
      if (!f) { if (saved.enabled !== false) missing++; continue; }
      for (const k of EDITABLE) if (saved[k] !== undefined) f[k] = saved[k] || f.base[k];
      f.enabled = saved.enabled !== false;
      order.push(f);
    }
    const used = new Set(order);
    S.files.forEach(f => { if (!used.has(f)) f.enabled = false; });
    S.files = [...order, ...S.files.filter(f => !used.has(f))];
  }
  renderAll(); changed();
  if (!quiet) toast(missing ? `Configuration loaded; ${missing} file(s) are not in this project.` : 'Configuration loaded');
  else if (missing) toast(`${missing} file(s) used by this figure are no longer in the project.`);
}

/* ------------------------------------------------------------------ project loading */
async function loadProject(id) {
  const info = await api(`/api/projects/${encodeURIComponent(id)}`);
  S.project = info.project;
  S.role = info.role;
  S.templates = info.templates;
  const saved = store.get(`${STORE}:${id}`, null);
  const overrides = saved?.files || {};
  let files = info.files.map(f => hydrate(f, overrides[f.id] || {}));
  if (saved?.order) {
    const rank = new Map(saved.order.map((fid, i) => [fid, i]));
    files = files.sort((a, b) => (rank.get(a.id) ?? 1e6) - (rank.get(b.id) ?? 1e6));
  }
  S.files = files;
  S.category = saved?.category ?? 1;
  S.plot = saved?.plot ?? 2;
  if (!S.meta.catalog[S.category]?.[1][S.plot]) { S.category = 1; S.plot = 2; }
  S.settings = saved?.settings || { ...(S.project.defaults || {}) };
  S.advanced = saved?.advanced || '';
  S.history = []; S.current = null; S.figureId = null; S.openFile = null;
  $('#figure').hidden = true;
  $('#interactive').hidden = true; $('#viewtoggle').hidden = true; $('#paper').classList.remove('interactive-on');
  const ph = $('#placeholder');
  ph.hidden = false; ph.classList.remove('error');
  $('#savebtn').disabled = true;
  renderAll(); renderHistory(); persist();
  await loadFigures();
  status(`${S.project.name}: ${S.files.length} files${readOnly() ? ' · view only' : ''}.`);
  if (!readOnly() && requirement().ok) run();
}

/* ------------------------------------------------------------------ wiring */
function renderAll() { renderTopbar(); renderProjectInfo(); renderCategories(); renderPlots(); renderFiles(); renderNeeds(); renderSettings(); renderCrumbs(); renderTemplates(); }
const MENUS = [['#exportbtn', '#exportmenu'], ['#configbtn', '#configmenu'], ['#tplbtn', '#tplmenu'], ['#userbtn', '#usermenu']];
function closeMenus(except) {
  for (const [btn, menu] of MENUS) {
    if (menu === except) continue;
    $(menu).hidden = true; $(btn).setAttribute('aria-expanded', 'false');
  }
}
function toggleMenu(btn, menu) {
  const open = $(menu).hidden;
  closeMenus(menu);
  $(menu).hidden = !open; $(btn).setAttribute('aria-expanded', String(open));
}
function setDataTab(tab) {
  $$('#datatabs button').forEach(b => b.setAttribute('aria-selected', String(b.dataset.tab === tab)));
  for (const t of ['files', 'figures']) $('#tab-' + t).hidden = t !== tab;
}
function wire() {
  $('#run').onclick = () => run();
  $('#savebtn').onclick = saveFigure;
  $('#live').onchange = e => { S.live = e.target.checked; persist(); if (S.live && S.resultRevision !== S.revision && requirement().ok) run(); };
  $('#project').onchange = e => { persist(); loadProject(e.target.value).catch(x => status(x.message, 'error')); };
  $('#theme').onclick = () => {
    const root = document.documentElement;
    const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    root.dataset.theme = dark ? 'light' : 'dark';
    store.set('gc-theme', root.dataset.theme);
  };
  $$('#datatabs button').forEach(b => { b.onclick = () => setDataTab(b.dataset.tab); });
  const drop = $('#dropzone'), input = $('#fileinput');
  drop.onclick = () => input.click();
  drop.onkeydown = e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); input.click(); } };
  input.onchange = () => { uploadAll(input.files); input.value = ''; };
  let depth = 0;
  const hasFiles = e => [...(e.dataTransfer?.types || [])].includes('Files');
  window.addEventListener('dragenter', e => { if (hasFiles(e) && !readOnly()) { depth++; document.body.classList.add('dragging'); } });
  window.addEventListener('dragleave', e => { if (hasFiles(e) && --depth <= 0) { depth = 0; document.body.classList.remove('dragging'); } });
  window.addEventListener('dragover', e => { if (hasFiles(e)) e.preventDefault(); });
  window.addEventListener('drop', e => {
    if (!hasFiles(e)) return;
    e.preventDefault(); depth = 0; document.body.classList.remove('dragging');
    setDataTab('files'); uploadAll(e.dataTransfer.files);
  });
  $$('#zoom button').forEach(b => { b.onclick = () => setZoom(b.dataset.zoom); });
  $$('#viewtoggle button').forEach(b => { b.onclick = () => { S.view = b.dataset.view; store.set('gc-view', S.view); setView(S.view); }; });
  $('#canvas').addEventListener('wheel', e => { if ((e.ctrlKey || e.metaKey) && S.current) { e.preventDefault(); setZoom(e.deltaY < 0 ? 'in' : 'out'); } }, { passive: false });
  $('#bgtoggle').onclick = e => { S.checker = !S.checker; $('#paper').classList.toggle('checker', S.checker); e.target.setAttribute('aria-pressed', String(S.checker)); };
  $('#exportbtn').onclick = e => { e.stopPropagation(); renderExportMenu(); toggleMenu('#exportbtn', '#exportmenu'); };
  $('#configbtn').onclick = e => { e.stopPropagation(); toggleMenu('#configbtn', '#configmenu'); };
  $('#tplbtn').onclick = e => { e.stopPropagation(); renderTemplates(); toggleMenu('#tplbtn', '#tplmenu'); };
  $('#userbtn').onclick = e => { e.stopPropagation(); toggleMenu('#userbtn', '#usermenu'); };
  $('#configmenu').onclick = async e => {
    const act = e.target.closest('button')?.dataset.act;
    if (!act) return;
    e.stopPropagation(); closeMenus();
    if (act === 'save') {
      const url = URL.createObjectURL(new Blob([JSON.stringify(configObject(), null, 2)], { type: 'application/json' }));
      h('a', { href: url, download: slug(val('title')) + '.config.json' }).click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    } else if (act === 'load') $('#configinput').click();
    else if (act === 'copy') { await navigator.clipboard.writeText(JSON.stringify(configObject(), null, 2)); toast('Configuration copied'); }
  };
  $('#configinput').onchange = async e => {
    const file = e.target.files[0];
    e.target.value = '';
    if (!file) return;
    try { applyConfig(JSON.parse(await file.text())); } catch (x) { toast('Not a valid configuration: ' + x.message); }
  };
  $('#resetsettings').onclick = () => {
    FIELDS.filter(applies).map(f => f.k).filter(k => !REGION_KEYS.includes(k)).forEach(k => delete S.settings[k]);
    for (const [k, v] of Object.entries(S.project.defaults || {})) if (!(k in S.settings)) S.settings[k] = v;
    renderSettings(); renderFiles(); renderCrumbs(); changed();
    toast('Settings for this plot reset to project defaults (region and data kept)');
  };
  $('#settingsearch').oninput = debounce(e => { S.search = e.target.value; renderSettings(); }, 120);
  document.addEventListener('click', e => { if (!e.target.closest('.menu-wrap')) closeMenus(); });
  document.addEventListener('keydown', e => {
    if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') { e.preventDefault(); document.activeElement?.dispatchEvent?.(new Event('change', { bubbles: true })); setTimeout(run, 0); }
    else if (e.key === 'Escape') closeMenus();
    else if (e.key === '/' && !/INPUT|TEXTAREA|SELECT/.test(document.activeElement?.tagName)) { e.preventDefault(); $('#settingsearch').focus(); }
  });
}

async function boot() {
  S.view = store.get('gc-view', 'interactive');
  const theme = store.get('gc-theme', null);
  if (theme) document.documentElement.dataset.theme = theme;
  wire();
  try {
    [S.me, S.meta] = await Promise.all([api('/api/me'), api('/api/meta')]);
  } catch (e) {
    status('Cannot reach the studio server: ' + e.message, 'error');
    fill($('#placeholder'), h('h3', {}, 'Server offline'), h('p', {}, 'Start it with gcanvas-studio serve (or python3 web/server.py), then reload.'));
    return;
  }
  $('#engine').textContent = `${S.meta.engine} Studio · ${S.me.user.email === 'local' ? 'local mode' : S.me.user.email}`;
  const global = store.get(STORE, {});
  S.live = !!global.live;
  S.sections = { ...S.sections, ...(global.sections || {}) };
  $('#live').checked = S.live;
  S.projects = S.me.projects;
  if (!S.projects.length) {
    renderTopbar();
    fill($('#placeholder'), h('h3', {}, 'No projects yet'), h('p', {}, 'The core facility has not shared a project with this account. Contact them to get access.'));
    status('No projects available for ' + S.me.user.email);
    return;
  }
  const wanted = new URLSearchParams(location.search).get('project') || global.project;
  const first = S.projects.find(p => p.id === wanted) || S.projects.find(p => !p.demo) || S.projects[0];
  S.project = first;
  await loadProject(first.id);
}
boot();
