"""Independent B0[eps] references: python -m pip install mpmath; run this file.

Integrate the Feynman parameter logarithm at 90 decimal digits, splitting at
real zeros. No dilogarithm formula or OneLOop implementation is used.
"""
from pathlib import Path
import mpmath as mp

mp.mp.dps = 90
POINTS = [
    ('-3', '2', '0', '3', '0', '1'),
    ('10', '2', '0', '3', '0', '1'),
    ('20', '2', '0', '3', '0', '2'),
    ('0', '2', '0', '3', '0', '1'),
    ('0', '2', '0', '2', '0', '1'),
    ('0', '0', '0', '2', '0', '1'),
    ('0', '0', '0', '0', '0', '1'),
    ('-3', '0', '0', '0', '0', '2'),
    ('3', '0', '0', '0', '0', '2'),
    ('-3', '0', '0', '2', '0', '1'),
    ('3', '0', '0', '2', '0', '1'),
    ('2', '0', '0', '2', '0', '1'),
    ('4', '1', '0', '1', '0', '1'),
    ('1', '1', '0', '4', '0', '1'),
    ('1', '4', '0', '1', '0', '1'),
    ('-3', '-2', '0', '-3', '0', '1'),
    ('3', '-2', '0', '-3', '0', '1'),
    ('-3', '2', '-0.1', '3', '-0.2', '1'),
    ('20', '2', '-0.1', '3', '-0.2', '2'),
    ('-3', '-2', '-0.1', '3', '-0.2', '1'),
    ('-3', '2', '-0.1', '-3', '-0.2', '1'),
    ('3', '0', '0', '2', '-0.1', '1'),
    ('0', '2', '-0.1', '3', '-0.2', '1'),
    ('0.01', '2', '0', '3', '0', '1'),
    ('0.01', '2', '-0.1', '3', '-0.2', '1'),
    ('-3', '2', '0', '2.0001', '0', '1'),
    ('3.99', '1', '0', '1', '0', '1'),
    ('4.01', '1', '0', '1', '0', '1'),
]

def reference(p, a, b, scale):
    if not any((p,a,b)):
        return mp.mpc(0)
    if p == 4 and a == b == scale == 1:
        # Q=(2x-1)^2: integral log^2(Q) = 8 exactly.
        return mp.zeta(2) + 4
    roots = mp.polyroots([p, b-a-p, a], maxsteps=1000) if p else [a/(a-b)] if a!=b else []
    cuts = sorted({mp.mpf(0),mp.mpf(1),*(mp.re(r)+t*abs(mp.im(r)) for r in roots for t in [-2,-1,0,1,2] if 0<mp.re(r)+t*abs(mp.im(r))<1)})
    def f(x):
        z=(a*(1-x)+b*x-p*x*(1-x))/scale
        if not z: return mp.mpf(0) # isolated endpoint; measure zero
        l=mp.log(z) if mp.im(z) or mp.re(z)>0 else mp.log(-z)-mp.pi*1j
        return l*l
    return mp.zeta(2)+mp.quad(f,cuts)/2

rows=['# p m0.re m0.im m1.re m1.im mu2 B0eps.re B0eps.im']
for row in POINTS:
    p, ar, ai, br, bi, scale=map(mp.mpf,row)
    value=reference(p,mp.mpc(ar,ai),mp.mpc(br,bi),scale)
    rows.append(' '.join((*row,mp.nstr(value.real,75),mp.nstr(value.imag,75))))
Path(__file__).with_name('b0_epsilon.tsv').write_text('\n'.join(rows)+'\n')
print(f'Generated {len(POINTS)} independent B0 epsilon fixtures')
