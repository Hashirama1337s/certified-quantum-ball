#!/usr/bin/env python
"""make_scorecard.py -- build scorecard.svg and SCORECARD.md: every sealed test of this work, with every sealed item,
failures included.

The five sealed bars of this work (QB1, QB1b, QB2, QB2b, QB2c) and the three inherited inputs (E5, M, BG2) are
listed with their outcomes; each outcome, and each per-item outcome where the ledger row records one, is checked against
the ledger export results/ledger-quantum.json (its settle rows). The script stops if anything disagrees.
Run from anywhere:  python tools/make_scorecard.py
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

H, F = 'HOLD', 'FAILED'
# (test, seal, what, outcome, source, [(item, outcome, what as scored)])
BARS = [
    ('QB1', '54caa4012c', 'Certified Bach source on Evans-Coleman (first instrument)', F, 'results/seals/BAR-QB1.md', [
        ('PQ1', H, 'sonic point not Bach-flat: b_uu = 0.44667377966953663 +/- 1.3e-17'),
        ('PQ2', F, 'covering sign map: signs certified at 691 grid points only, not on covering cells'),
        ('PQ3', H, 'centre law: c0 = -173.10044439229957 +/- 3e-13'),
        ('PQ4', H, 'routes A and B: 28 comparisons at 19 distinct x overlap'),
        ('PQ5', H, 'controls: Hayward, Reissner-Nordstrom, Minkowski, k = 1/2 Friedmann, sign-flip null'),
        ('PQ6', F, 'off-solution null ill-posed: the shifted background still solves the ODE locally'),
        ('PQ7', F, 'D_new(2) certified at x_m = -4.6, not at the sealed x_m = -2')]),
    ('QB1b', '4d3b12ef39', 'Re-test: Taylor-model sign map, off-shell null, bridge, seam', H, 'results/seals/BAR-QB1b.md', [
        ('PQ2b', H, '2261 cells cover [-12, 12]; exactly one zero, width 3.7e-9'),
        ('PQ6b', H, 'omega x (1 + 1e-12): residual excludes 0 at 24 of 28 points'),
        ('PQ7b', H, 'D_new(2) = -2.4126270857981 +/- 1.1e-9 at x_m = -4.6'),
        ('PQ8', H, 'centre series and tube overlap at x = -9/2')]),
    ('QB2', 'f2bcfadb3e', 'First-order response at rate 2 and Psi2', F, 'results/seals/BAR-QB2.md', [
        ('PQ2-1', H, 'source Taylor models match QB1b'),
        ('PQ2-2', H, 'four HKA source rows exact'),
        ('PQ2-3', F, 'conservation: e encloses 0 but width 2e-5 > 1e-10; (0, 12] not covered'),
        ('PQ2-4', F, 'centre coefficients <= 3.5e-16, not the sealed 1e-30'),
        ('PQ2-5', F, 'controls partial: C1 at x = 0.5, 3 not run; C8 theta-theta residual not computed'),
        ('PQ2-6', F, 'nulls partial: N-a refuses, but the design formula for |e| was wrong (1e-4 rel)'),
        ('PQ2-7', H, 'x_V unique in U, V\' > 0'),
        ('PQ2-8', H, 'headline: Psi2 = 3.4616358037141 +/- 1.15e-7 (ledger wording), positive')]),
    ('QB2b', '373f306e9f', 'Closing the QB2 gaps (sharper BG2 inward tube)', F, 'results/seals/BAR-QB2b.md', [
        ('PQ2b-1', F, 'e encloses 0, width 6.66e-10 > 1e-10 (D_t <= 3.4e-13 passes)'),
        ('PQ2b-2', H, 'exact centre orders; c0 = -(2/3)W0(8W0-3)e^{-4xc}'),
        ('PQ2b-3', H, 'domain lemma (x <= 0 only); C1 at five points'),
        ('PQ2b-4', H, 'C8 theta-theta, C8\', N-a (corrected), N-c, C2 as predicted'),
        ('PQ2b-5', H, 'headline replay: Psi2 = 3.4616358037141657 +/- 3.0e-12 (ledger wording; certified ball ± 3.02e-12); route B overlaps')]),
    ('QB2c', '7aab7efe5a', 'Conservation width, theta-theta over the domain, replication row', F, 'results/seals/BAR-QB2c.md', [
        ('PQ2c-1', F, 'not run as sealed: a probe showed the tight stage cannot move the width; the QB2b value '
                      '(e encloses 0, width 6.66e-10) is replayed'),
        ('PQ2c-2', H, 'theta-theta residual encloses 0 at five points; N-b refuses at four, and at x = -4.6 only on '
                      'the centre-side representation chosen POST HOC (the sealed sonic-side one does not exclude 0)'),
        ('PQ2c-3', H, 'headline replay: Psi2 = 3.4616358037141657 +/- 3.0e-12 (ledger wording; certified ball ± 3.02e-12)'),
        ('PQ2c-4', H, 'blind replication converts as 30 pi x_s^2 Psi_blind / A(x_V); rel <= 6.9e-12 (POST HOC conversion)')]),
]
INPUTS = [
    ('E5', 'd0af75803f', 'Existence of the Evans-Coleman solution (published, DOI 10.5281/zenodo.24303877)', H),
    ('M', 'b1b0ee1cde', 'Exactly one growing spherical mode (published, same record)', H),
    ('BG2', 'ed61d09268', 'Sharper sonic box and outward record (not in any published record; shipped here)', H),
]

FONT = 'Segoe UI, Helvetica, Arial, sans-serif'
INK, MUTED, ACCENT, BAND, LINE, BORDER = '#e8edf2', '#9aa4ae', '#39ff14', '#262b31', '#22272e', '#2f353d'
CHIP = {'HOLD': ('#12301c', '#39ff14'), 'FAILED': ('#3a1a1a', '#ff6b6b'), 'VOID': ('#33290f', '#ffc94d')}


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def check():
    rows = json.load(open(os.path.join(ROOT, 'results', 'ledger-quantum.json'), encoding='utf-8'))
    st = {r['seal_id']: r for r in rows if r.get('kind') == 'settle'}
    for tid, seal, _, out, _, items in BARS:
        if st.get(seal, {}).get('outcome') != out:
            sys.exit('ledger export disagrees for %s %s' % (tid, seal))
        led = st[seal].get('items')
        if led:
            mine = {i: o for i, o, _ in items}
            if led != mine:
                sys.exit('ledger items disagree for %s: %r vs %r' % (tid, led, mine))
    for tid, seal, _, out in INPUTS:
        if st.get(seal, {}).get('outcome') != out:
            sys.exit('ledger export disagrees for input %s %s' % (tid, seal))


def chip(x, y, outcome, small=False):
    fill, fg = CHIP[outcome]
    w = (len(outcome) * (6.6 if small else 7.4)) + (14 if small else 20)
    h = 16 if small else 20
    fs = 10 if small else 11.5
    return (f'<rect x="{x}" y="{y}" width="{w:.0f}" height="{h}" rx="{h / 2:.0f}" fill="{fill}"/>'
            f'<text x="{x + (7 if small else 10)}" y="{y + h - (4 if small else 6)}" font-size="{fs}" font-weight="600" '
            f'fill="{fg}">{outcome}</text>'), w


def build():
    check()
    nh = sum(1 for b in BARS if b[3] == H)
    items = [i for b in BARS for i in b[5]]
    ih = sum(1 for i in items if i[1] == H)
    title = 'Scorecard: %d sealed bars (%d HOLD, %d FAILED by item); %d sealed items (%d HOLD, %d FAILED)' % (
        len(BARS), nh, len(BARS) - nh, len(items), ih, len(items) - ih)
    sub = 'Sealed = predictions and pass bars hashed before the run. Failures are shown as failures. Details: SCORECARD.md.'
    W, y, body = 980, 66, []
    body.append(f'<rect x="10" y="{y}" width="{W - 20}" height="22" fill="{BAND}"/>')
    body.append(f'<text x="20" y="{y + 16}" font-size="11.5" font-weight="600" fill="{ACCENT}" letter-spacing="0.5">'
                f'THIS WORK: THE QUANTUM SOURCE (QB1) AND THE RESPONSE (QB2)</text>')
    y += 22
    for tid, seal, what, out, _, its in BARS:
        body.append(f'<text x="20" y="{y + 22}" font-size="12.5" font-family="Consolas, Menlo, monospace" fill="{INK}">'
                    f'{esc(tid + " " + seal)}</text>')
        body.append(f'<text x="190" y="{y + 22}" font-size="12.5" fill="{INK}">{esc(what)}</text>')
        c, _ = chip(W - 120, y + 8, out)
        body.append(c)
        y += 30
        x = 190
        for it, o, _ in its:
            if x + len(it) * 6.6 + 80 > W - 20:
                x, y = 190, y + 22
            body.append(f'<text x="{x}" y="{y + 12}" font-size="10.5" font-family="Consolas, Menlo, monospace" '
                        f'fill="{MUTED}">{esc(it)}</text>')
            c, w = chip(x + len(it) * 6.6 + 6, y, o, small=True)
            body.append(c)
            x += len(it) * 6.6 + 6 + w + 14
        y += 24
        body.append(f'<line x1="10" y1="{y + 1}" x2="{W - 10}" y2="{y + 1}" stroke="{LINE}"/>')
    y += 8
    body.append(f'<rect x="10" y="{y}" width="{W - 20}" height="22" fill="{BAND}"/>')
    body.append(f'<text x="20" y="{y + 16}" font-size="11.5" font-weight="600" fill="{ACCENT}" letter-spacing="0.5">'
                f'INHERITED INPUTS (EARLIER SEALED WORK)</text>')
    y += 22
    for tid, seal, what, out in INPUTS:
        body.append(f'<text x="20" y="{y + 22}" font-size="12.5" font-family="Consolas, Menlo, monospace" fill="{INK}">'
                    f'{esc(tid + " " + seal)}</text>')
        body.append(f'<text x="190" y="{y + 22}" font-size="12.5" fill="{INK}">{esc(what)}</text>')
        c, _ = chip(W - 120, y + 8, out)
        body.append(c)
        y += 30
        body.append(f'<line x1="10" y1="{y + 1}" x2="{W - 10}" y2="{y + 1}" stroke="{LINE}"/>')
    Ht = y + 14
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Ht}" viewBox="0 0 {W} {Ht}" font-family="{FONT}" '
           f'role="img" aria-label="{esc(title)}">\n<title>{esc(title)}</title>\n'
           f'<rect width="{W}" height="{Ht}" rx="10" fill="#1c2026" stroke="{BORDER}"/>\n'
           f'<text x="20" y="34" font-size="17" font-weight="600" fill="{INK}">{esc(title)}</text>\n'
           f'<text x="20" y="54" font-size="12" fill="{MUTED}">{esc(sub)}</text>\n' + '\n'.join(body) + '\n</svg>\n')
    with open(os.path.join(ROOT, 'scorecard.svg'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(svg)
    ch = {k: '![%s](brand/chips/%s.svg)' % (k, k.lower()) for k in CHIP}
    md = ['# Scorecard', '',
          '**In plain words.** Before each test ran, we wrote down what it had to show and hashed that text, so it could not',
          'be changed afterwards. Each test is scored exactly as written: one failed item makes the whole test FAILED. This',
          'page lists every such test of this work and every item in it, the failures included. The headline number held',
          'in all three tests that checked it; the conservation-width bar never held.', '',
          title + '.', '',
          'Generated by `tools/make_scorecard.py`. A sealed test is one whose predictions and pass bars were hashed before',
          'its run, and it is scored as written: a bar HOLDs only if every item holds. Every outcome below, and every',
          'per-item outcome the ledger records, is checked against the ledger export `results/ledger-quantum.json`.', '',
          '| Test | Seal | What was predicted | Verdict | Source |', '|---|---|---|---|---|']
    for tid, seal, what, out, src, _ in BARS:
        md.append(f'| {tid} | `{seal}` | {what} | {ch[out]} | [{src}]({src}) |')
    md += ['', '## Every sealed item', '', '| Test | Item | Verdict | As scored |', '|---|---|---|---|']
    for tid, _, _, _, _, its in BARS:
        for it, o, d in its:
            md.append(f'| {tid} | {it} | {ch[o]} | {d} |')
    md += ['', '## Inherited inputs', '', '| Test | Seal | What | Verdict |', '|---|---|---|---|']
    for tid, seal, what, out in INPUTS:
        md.append(f'| {tid} | `{seal}` | {what} | {ch[out]} |')
    md += ['', 'Notes:',
           '- The headline Psi2 held in all three QB2 bars (PQ2-8, PQ2b-5, PQ2c-3). The sealed conservation-width bar',
           '  (e width <= 1e-10) FAILED in QB2 and QB2b (e encloses 0, width 6.66e-10 in QB2b). In QB2c it was not run as',
           '  sealed: a probe showed the tight stage cannot move the width, and the QB2b value is replayed (FAILED). Our',
           '  estimate of the floor, from that probe and not a certificate: the step balls of the pinned inward background',
           '  tube (Taylor order KT = 36) times the state sensitivity of the transported system (`code/quantum/seats/build-5/RUN.md`).',
           '- The QB1b settle row carries an automatic pin warning ("STALE: result_qb1b.json names instrument hash',
           '  582e24236c8f..."). That hash is `engine/cert/lin_stages.py`, recorded by the instrument as provenance; the file',
           '  is unchanged and pinned under b1b0ee1cde. It is not a changed instrument.',
           '- QB1 PQ6 was an ill-posed null (an off-solution sonic velocity still solves the ODE locally). QB1b PQ6b replaced it',
           '  by a real off-shell null, which refuses.',
           '- An amendment sealed under QB2: 4cf7213207 (`results/seals/AMENDMENT-QB2-v3.md`, text fixes only).',
           '- The blind replication row (PQ2c-4) compares a double-precision blind computation with no enclosure; it is a',
           '  replication, not a certificate. The conversion factor 30 pi x_s^2 / A(x_V) was found after both numbers existed.']
    text = '\n'.join(md) + '\n'
    sys.path.insert(0, os.path.join(ROOT, 'tools'))
    import fill_balls                      # every printed ball verbatim or rounded outward from its stored string
    text, rows = fill_balls.fix_text(text)
    with open(os.path.join(ROOT, 'SCORECARD.md'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(text)
    print(title)


if __name__ == '__main__':
    build()
