'use strict';
/* GenomeCanvas Studio admin console: projects, membership, users, templates, activity. */
const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];
function h(tag, attrs = {}, ...kids) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v == null || v === false) continue;
    if (k === 'class') el.className = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else if (k === 'value') el.value = v;
    else if (k === 'checked') el.checked = !!v;
    else el.setAttribute(k, v === true ? '' : v);
  }
  for (const kid of kids.flat()) if (kid != null && kid !== false) el.append(kid instanceof Node ? kid : String(kid));
  return el;
}
const fmtBytes = n => !n ? '0' : n < 1048576 ? (n / 1024).toFixed(0) + ' KB' : n < 1073741824 ? (n / 1048576).toFixed(1) + ' MB' : (n / 1073741824).toFixed(2) + ' GB';
const fmtDate = t => t ? new Date(t * 1000).toLocaleString(undefined, { dateStyle: 'medium', timeStyle: 'short' }) : '—';
let toastTimer;
function toast(msg) { const t = $('#toast'); t.textContent = msg; t.hidden = false; clearTimeout(toastTimer); toastTimer = setTimeout(() => { t.hidden = true; }, 4000); }

async function api(path, opt = {}) {
  const res = await fetch(path, { ...opt, credentials: 'same-origin', headers: { 'X-GenomeCanvas': 'studio', ...(opt.body ? { 'Content-Type': 'application/json' } : {}) } });
  if (res.status === 401) { location.href = '/login'; throw new Error('Please sign in'); }
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || 'Request failed');
  return data;
}
const send = (path, method, body) => api(path, { method, body: JSON.stringify(body || {}) });

/** Modal form; resolves to values or null. */
function ask(title, fields, { ok = 'Save', extra } = {}) {
  return new Promise(resolve => {
    const dlg = $('#dialog');
    $('#dialogtitle').textContent = title;
    $('#dialogok').textContent = ok;
    $('#dialogok').hidden = !fields.length;
    const inputs = {};
    $('#dialogbody').replaceChildren(...(extra ? [extra] : []), ...fields.map(f => {
      let input;
      if (f.type === 'select') input = h('select', {}, ...f.options.map(([v, t]) => h('option', { value: v, selected: v === f.value }, t)));
      else if (f.type === 'textarea') input = h('textarea', { rows: 3, value: f.value || '' });
      else if (f.type === 'check') input = h('input', { type: 'checkbox', checked: !!f.value });
      else input = h('input', { type: f.type || 'text', value: f.value ?? '', placeholder: f.placeholder || '' });
      inputs[f.name] = input;
      return f.type === 'check' ? h('label', { class: 'dialog-check' }, input, f.label)
        : h('label', { class: 'dialog-field' }, h('span', {}, f.label), input, f.help ? h('small', {}, f.help) : null);
    }));
    dlg.onclose = () => resolve(dlg.returnValue === 'ok'
      ? Object.fromEntries(Object.entries(inputs).map(([k, el]) => [k, el.type === 'checkbox' ? el.checked : el.value.trim()])) : null);
    dlg.returnValue = '';
    dlg.showModal();
  });
}
function showInvite(result, email) {
  if (!result.invite_url) { if (result.note) toast(result.note); return; }
  const input = h('input', { class: 'input', value: result.invite_url, readonly: true });
  ask(`Invitation for ${email}`, [], {
    extra: h('div', { class: 'invite' },
      h('p', {}, 'Send this link to the user. It sets their password and is valid for 7 days.'), input,
      h('button', { class: 'btn small', type: 'button', onclick: async () => { await navigator.clipboard.writeText(result.invite_url); toast('Link copied'); } }, 'Copy link')),
  });
  setTimeout(() => input.select(), 50);
}

let data = null;
async function load() {
  data = await api('/api/admin/overview');
  $('#roots').textContent = `Data root: ${data.data_root || '(not set)'} · Shared annotations: ${data.shared_root || '(not set)'} · Sign-in: ${data.auth}`;
  renderProjects(); renderUsers(); renderTemplates();
}

function renderProjects() {
  const q = $('#projectfilter').value.trim().toLowerCase();
  const items = data.projects.filter(p => !q || [p.id, p.name, ...p.members.map(m => m.email)].join(' ').toLowerCase().includes(q));
  const today = new Date().toISOString().slice(0, 10);
  $('#projects').replaceChildren(...items.map(p => {
    const expired = p.expires && p.expires < today;
    const state = p.status !== 'active' ? 'archived' : expired ? 'expired' : 'active';
    return h('article', { class: 'admin-card' + (state !== 'active' ? ' inactive' : '') },
      h('header', {},
        h('div', {}, h('h3', {}, p.name), h('code', {}, p.id), ' ', h('span', { class: 'chip ' + (state === 'active' ? 'ok' : 'miss') }, state),
          p.genome ? h('span', { class: 'chip' }, p.genome) : null, p.demo ? h('span', { class: 'chip' }, 'demo · open to all users') : null),
        h('a', { class: 'btn small ghost', href: `/?project=${encodeURIComponent(p.id)}` }, 'Open')),
      p.description ? h('p', { class: 'hint' }, p.description) : null,
      h('dl', { class: 'stats' },
        h('div', {}, h('dt', {}, 'Delivered files'), h('dd', {}, p.files?.delivered || 0)),
        h('div', {}, h('dt', {}, 'Uploads'), h('dd', {}, `${p.files?.uploads || 0} · ${fmtBytes(p.upload_bytes)}${p.quota_bytes ? ' / ' + fmtBytes(p.quota_bytes) : ''}`)),
        h('div', {}, h('dt', {}, 'Figures (30 d)'), h('dd', {}, `${p.usage.jobs || 0}${p.usage.errors ? ` · ${p.usage.errors} failed` : ''}`)),
        h('div', {}, h('dt', {}, 'Saved figures'), h('dd', {}, p.figures)),
        h('div', {}, h('dt', {}, 'Folder'), h('dd', {}, p.data_dir || '—')),
        h('div', {}, h('dt', {}, 'Access until'), h('dd', {}, p.expires || 'no end date'))),
      p.demo ? null : h('div', { class: 'members' },
        h('strong', {}, `Members (${p.members.length})`),
        ...p.members.map(m => h('span', { class: 'member' }, m.email, h('em', {}, m.role),
          h('button', { title: 'Remove from project', onclick: () => setMember(p, m.email, null) }, '×'))),
        h('button', { class: 'btn small', onclick: () => addMember(p) }, '+ Add member')),
      h('footer', { class: 'inline-btns' },
        p.data_dir ? h('button', { class: 'btn small', onclick: () => rescan(p) }, 'Rescan data') : null,
        p.demo ? null : h('button', { class: 'btn small ghost', onclick: () => editProject(p) }, 'Edit'),
        p.demo ? null : h('button', { class: 'btn small ghost', onclick: () => send(`/api/admin/projects/${p.id}`, 'PATCH', { status: p.status === 'active' ? 'archived' : 'active' }).then(load).catch(e => toast(e.message)) }, p.status === 'active' ? 'Archive' : 'Re-activate')));
  }));
  if (!items.length) $('#projects').append(h('p', { class: 'empty' }, 'No projects match.'));
}
const projectFields = p => [
  p ? null : { name: 'id', label: 'Project ID', placeholder: 'e.g. P2026-017-Smith', help: 'Letters, digits, . _ - ; usually your LIMS / order number.' },
  { name: 'name', label: 'Display name', value: p?.name, placeholder: 'Smith lab · H3K27ac ChIP-seq' },
  { name: 'data_dir', label: 'Data folder (relative to data root)', value: p?.data_dir, placeholder: 'smith/P2026-017', help: 'Put samples.tsv (+ optional project.yaml) here to control labels, groups and colours.' },
  { name: 'genome', label: 'Genome', value: p?.genome, placeholder: 'hg38' },
  { name: 'description', label: 'Description', type: 'textarea', value: p?.description },
  { name: 'expires', label: 'Access until', type: 'date', value: p?.expires || '', help: 'Leave empty for no end date.' },
  { name: 'quota_gb', label: 'Upload quota (GB)', type: 'number', value: p ? (p.quota_bytes / 1024 ** 3).toFixed(1) : '' },
].filter(Boolean);
async function newProject() {
  const v = await ask('New project', projectFields(null), { ok: 'Create' });
  if (!v) return;
  try {
    const r = await send('/api/admin/projects', 'POST', v);
    toast(`Created ${r.project.id}${r.scan?.files != null ? ` · ${r.scan.files} files found` : ''}${r.scan?.warnings?.length ? ' · ' + r.scan.warnings.join('; ') : ''}`);
    await load();
  } catch (e) { toast(e.message); }
}
async function editProject(p) {
  const v = await ask(`Edit ${p.id}`, projectFields(p));
  if (!v) return;
  try { await send(`/api/admin/projects/${p.id}`, 'PATCH', v); await load(); toast('Saved'); } catch (e) { toast(e.message); }
}
async function rescan(p) {
  try {
    const { scan } = await send(`/api/admin/projects/${p.id}/scan`, 'POST');
    toast(`${p.id}: ${scan.files ?? 0} files (${scan.added || 0} new, ${scan.removed || 0} removed)${scan.warnings?.length ? ' · ' + scan.warnings.join('; ') : ''}`);
    await load();
  } catch (e) { toast(e.message); }
}
async function addMember(p) {
  const v = await ask(`Add member to ${p.id}`, [
    { name: 'email', label: 'E-mail', type: 'email', help: 'New addresses get an account and an invitation link.' },
    { name: 'role', label: 'Access', type: 'select', value: 'editor', options: [['editor', 'Editor: plot, upload, save figures'], ['viewer', 'Viewer: open and download saved figures']] },
  ], { ok: 'Add' });
  if (v?.email) setMember(p, v.email, v.role);
}
async function setMember(p, email, role) {
  try {
    const r = await send(`/api/admin/projects/${p.id}/members`, 'POST', { email, role });
    await load();
    if (r.invite_url) showInvite(r, email); else toast(role ? `${email} is now ${role} on ${p.id}` : `Removed ${email}`);
  } catch (e) { toast(e.message); }
}

function renderUsers() {
  const q = $('#userfilter').value.trim().toLowerCase();
  const memberships = {};
  for (const p of data.projects) for (const m of p.members) (memberships[m.id] ||= []).push(`${p.id} (${m.role})`);
  const rows = data.users.filter(u => !q || `${u.email} ${u.name}`.toLowerCase().includes(q));
  $('#users').replaceChildren(
    h('thead', {}, h('tr', {}, ...['User', 'Projects', 'Last sign-in', 'Status', ''].map(t => h('th', {}, t)))),
    h('tbody', {}, ...rows.map(u => h('tr', { class: u.active ? '' : 'inactive' },
      h('td', {}, h('strong', {}, u.name || ''), h('div', { class: 'hint' }, u.email)),
      h('td', {}, (memberships[u.id] || []).join(', ') || (u.is_admin ? 'all (admin)' : '—')),
      h('td', {}, fmtDate(u.last_login)),
      h('td', {}, u.is_admin ? h('span', { class: 'chip ok' }, 'admin') : null, u.active ? null : h('span', { class: 'chip miss' }, 'disabled')),
      h('td', { class: 'actions' },
        data.auth === 'local' ? h('button', { class: 'btn small ghost', onclick: async () => showInvite(await send(`/api/admin/users/${u.id}/invite`, 'POST'), u.email) }, 'Reset link') : null,
        h('button', { class: 'btn small ghost', onclick: () => send(`/api/admin/users/${u.id}`, 'PATCH', { is_admin: !u.is_admin }).then(load).catch(e => toast(e.message)) }, u.is_admin ? 'Remove admin' : 'Make admin'),
        h('button', { class: 'btn small ghost', onclick: () => send(`/api/admin/users/${u.id}`, 'PATCH', { active: !u.active }).then(load).catch(e => toast(e.message)) }, u.active ? 'Disable' : 'Enable'))))));
}
async function newUser() {
  const v = await ask('New user', [
    { name: 'email', label: 'E-mail', type: 'email' }, { name: 'name', label: 'Name' },
    { name: 'is_admin', label: 'Core staff (administrator)', type: 'check' },
  ], { ok: 'Create' });
  if (!v) return;
  try { const r = await send('/api/admin/users', 'POST', v); await load(); showInvite(r, v.email); } catch (e) { toast(e.message); }
}

function renderTemplates() {
  $('#templates').replaceChildren(
    h('thead', {}, h('tr', {}, ...['Template', 'Scope', 'Plot', ''].map(t => h('th', {}, t)))),
    h('tbody', {}, ...data.templates.map(t => h('tr', {},
      h('td', {}, h('strong', {}, t.name), t.description ? h('div', { class: 'hint' }, t.description) : null),
      h('td', {}, t.scope === 'global' ? 'all projects' : 'project', t.managed ? h('span', { class: 'chip' }, 'project.yaml') : null),
      h('td', {}, `${t.config.category ?? ''}/${t.config.plot ?? ''}`),
      h('td', { class: 'actions' }, t.managed ? null : h('button', { class: 'btn small ghost danger', onclick: () => send(`/api/templates/${t.id}`, 'DELETE').then(load).catch(e => toast(e.message)) }, 'Delete'))))));
  if (!data.templates.length) $('#templates').append(h('tbody', {}, h('tr', {}, h('td', { colspan: 4, class: 'empty' }, 'No templates yet.'))));
}

let audit = [];
async function loadAudit() {
  audit = (await api('/api/admin/audit?limit=1000')).audit;
  renderAudit();
}
function renderAudit() {
  const q = $('#auditfilter').value.trim().toLowerCase();
  const rows = audit.filter(a => !q || `${a.email} ${a.project_id} ${a.action} ${a.detail}`.toLowerCase().includes(q)).slice(0, 500);
  $('#audit').replaceChildren(
    h('thead', {}, h('tr', {}, ...['Time', 'User', 'Project', 'Action', 'Detail', 'IP'].map(t => h('th', {}, t)))),
    h('tbody', {}, ...rows.map(a => h('tr', {}, h('td', {}, fmtDate(a.ts)), h('td', {}, a.email || '—'), h('td', {}, a.project_id || ''),
      h('td', {}, h('code', {}, a.action)), h('td', { class: 'detail' }, a.detail || ''), h('td', { class: 'hint' }, a.ip || '')))));
}

function setTab(tab) {
  $$('#tabs button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.tab === tab)));
  for (const t of ['projects', 'users', 'templates', 'activity']) $('#tab-' + t).hidden = t !== tab;
  if (tab === 'activity') loadAudit().catch(e => toast(e.message));
  history.replaceState(null, '', '#' + tab);
}
$$('#tabs button').forEach(b => { b.onclick = () => setTab(b.dataset.tab); });
$('#newproject').onclick = newProject;
$('#newuser').onclick = newUser;
$('#projectfilter').oninput = renderProjects;
$('#userfilter').oninput = renderUsers;
$('#auditfilter').oninput = renderAudit;
$('#cleanup').onclick = async () => { try { const { cleanup } = await send('/api/admin/cleanup', 'POST'); toast(`Removed ${cleanup.jobs} old results, ${cleanup.uploads} expired uploads, ${cleanup.cache} cache entries`); } catch (e) { toast(e.message); } };
$('#rescanshared').onclick = async () => { try { const { scan } = await send('/api/admin/rescan-shared', 'POST'); toast(`${scan.files} shared annotation files`); } catch (e) { toast(e.message); } };
try { const t = localStorage.getItem('gc-theme'); if (t) document.documentElement.dataset.theme = JSON.parse(t); } catch { /* ignore */ }
window.addEventListener('hashchange', () => setTab(location.hash.slice(1) || 'projects'));
load().then(() => setTab(location.hash.slice(1) || 'projects')).catch(e => toast(e.message));
