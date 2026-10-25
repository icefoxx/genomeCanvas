'use strict';
/* Sign-in and invitation acceptance for GenomeCanvas Studio (password accounts). */
const $ = s => document.querySelector(s);
const token = new URLSearchParams(location.search).get('token');
const inviting = location.pathname === '/invite' && token;

if (inviting) {
  $('#heading').textContent = 'Set your password';
  $('#lead').textContent = 'Welcome! Choose a password (at least 10 characters) to activate your account.';
  $('#emailrow').hidden = true;
  $('#confirmrow').hidden = false;
  $('#pwlabel').textContent = 'New password';
  $('#password').autocomplete = 'new-password';
  $('#submit').textContent = 'Activate account';
  $('#foot').hidden = true;
}

function fail(message) {
  const el = $('#error');
  el.textContent = message;
  el.hidden = false;
}

$('#form').addEventListener('submit', async event => {
  event.preventDefault();
  $('#error').hidden = true;
  const password = $('#password').value;
  let path, body;
  if (inviting) {
    if (password.length < 10) return fail('Use at least 10 characters.');
    if (password !== $('#confirm').value) return fail('The two passwords do not match.');
    path = '/api/invite';
    body = { token, password };
  } else {
    path = '/api/login';
    body = { email: $('#email').value.trim(), password };
  }
  $('#submit').disabled = true;
  try {
    const res = await fetch(path, {
      method: 'POST', credentials: 'same-origin',
      headers: { 'Content-Type': 'application/json', 'X-GenomeCanvas': 'studio' },
      body: JSON.stringify(body),
    });
    const data = await res.json().catch(() => ({}));
    if (!res.ok) throw new Error(data.error || 'Sign-in failed');
    location.href = '/';
  } catch (e) {
    fail(e.message);
  } finally {
    $('#submit').disabled = false;
  }
});
