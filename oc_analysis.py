import csv, math
from datetime import datetime, timedelta
from collections import defaultdict

DATA = '/Users/mgarrod/CascadeProjects/tcrb/TCrB_V_1946-2026.csv'
P = 227.5528
T0 = 2455828.9
HALF_WINDOW = 50.0   # days, < P/4=56.89, keeps clear of secondary minimum
MIN_NIGHTS = 10

def jd_to_decyear(jd):
    d = datetime(1970,1,1) + timedelta(days=jd-2440587.5)
    year_start = datetime(d.year,1,1)
    year_end = datetime(d.year+1,1,1)
    frac = (d - year_start).total_seconds() / (year_end - year_start).total_seconds()
    return d.year + frac

# ---- load & filter ----
rows = []
with open(DATA, newline='') as f:
    r = csv.DictReader(f)
    for row in r:
        if row['fainterthan'] == 'True':
            continue
        try:
            jd = float(row['jd'])
            mag = float(row['mag'])
        except ValueError:
            continue
        rows.append((jd, mag))

rows.sort()
print("Usable rows (fainterthan excluded):", len(rows))

# ---- nightly bin (per integer JD day) ----
nightly = defaultdict(list)
for jd, mag in rows:
    key = round(jd)  # night bucket
    nightly[key].append((jd, mag))

nights = []
for key, obs in nightly.items():
    jds = [o[0] for o in obs]
    mags = [o[1] for o in obs]
    nights.append((sum(jds)/len(jds), sum(mags)/len(mags), len(mags)))
nights.sort()
print("Nightly-binned points:", len(nights))

jd_min = nights[0][0]
jd_max = nights[-1][0]
E_min = math.floor((jd_min - T0)/P) - 1
E_max = math.ceil((jd_max - T0)/P) + 1
print("Epoch range to scan:", E_min, "to", E_max)

def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
           -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
           +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))

def inv3(m):
    D = det3(m)
    if abs(D) < 1e-12:
        return None
    cof = [[0]*3 for _ in range(3)]
    def minor(mm, i, j):
        rows_ = [r for k,r in enumerate(mm) if k!=i]
        return [[v for k,v in enumerate(r) if k!=j] for r in rows_]
    def det2(m2):
        return m2[0][0]*m2[1][1]-m2[0][1]*m2[1][0]
    for i in range(3):
        for j in range(3):
            cof[i][j] = ((-1)**(i+j)) * det2(minor(m,i,j))
    # inverse = adjugate/D, adjugate = transpose of cofactor matrix
    inv = [[cof[j][i]/D for j in range(3)] for i in range(3)]
    return inv

def weighted_quadratic_fit(t, y, w):
    # solve weighted least squares for y = a*t^2 + b*t + c  (t already centered)
    Sw = sum(w)
    Swt = sum(wi*ti for wi,ti in zip(w,t))
    Swt2 = sum(wi*ti*ti for wi,ti in zip(w,t))
    Swt3 = sum(wi*ti**3 for wi,ti in zip(w,t))
    Swt4 = sum(wi*ti**4 for wi,ti in zip(w,t))
    Swy = sum(wi*yi for wi,yi in zip(w,y))
    Swty = sum(wi*ti*yi for wi,ti,yi in zip(w,t,y))
    Swt2y = sum(wi*ti*ti*yi for wi,ti,yi in zip(w,t,y))
    # normal equations matrix for [a,b,c]
    A = [[Swt4, Swt3, Swt2],
         [Swt3, Swt2, Swt],
         [Swt2, Swt,  Sw ]]
    B = [Swt2y, Swty, Swy]
    D = det3(A)
    if abs(D) < 1e-12:
        return None
    sols = []
    for col in range(3):
        Ai = [row[:] for row in A]
        for r_ in range(3):
            Ai[r_][col] = B[r_]
        sols.append(det3(Ai)/D)
    a,b,c = sols

    # covariance: sigma^2 * (X^T W X)^-1, sigma^2 = weighted SSR / (n-3)
    Ainv = inv3(A)
    resid = [yi-(a*ti*ti+b*ti+c) for ti,yi in zip(t,y)]
    n = len(t)
    if n <= 3 or Ainv is None:
        return a,b,c,None
    ssr = sum(wi*ri*ri for wi,ri in zip(w,resid))
    sigma2 = ssr/(n-3)
    # Ainv rows/cols correspond to [a,b,c]
    var_a, cov_ab = Ainv[0][0]*sigma2, Ainv[0][1]*sigma2
    var_b = Ainv[1][1]*sigma2
    cov = {'var_a': var_a, 'var_b': var_b, 'cov_ab': cov_ab}
    return a,b,c,cov

results = []
for E in range(E_min, E_max+1):
    t_calc = T0 + E*P
    window = [(jd, mag) for jd, mag, n in nights if abs(jd - t_calc) <= HALF_WINDOW]
    if len(window) < MIN_NIGHTS:
        continue
    # need coverage on both sides
    before = [w for w in window if w[0] < t_calc]
    after  = [w for w in window if w[0] >= t_calc]
    if len(before) < 3 or len(after) < 3:
        continue

    # iterative sigma-clip quadratic fit (2 passes)
    pts = window[:]
    fit = None
    for it in range(3):
        t = [p[0]-t_calc for p in pts]
        y = [p[1] for p in pts]
        w = [1.0]*len(pts)
        fit = weighted_quadratic_fit(t, y, w)
        if fit is None:
            break
        a,b,c,cov = fit
        resid = [yi - (a*ti*ti+b*ti+c) for ti,yi in zip(t,y)]
        sd = (sum(r*r for r in resid)/len(resid))**0.5
        if sd == 0:
            break
        new_pts = [p for p,r in zip(pts,resid) if abs(r) <= 3*sd]
        if len(new_pts) == len(pts) or len(new_pts) < MIN_NIGHTS:
            pts = new_pts if len(new_pts) >= MIN_NIGHTS else pts
            break
        pts = new_pts

    if fit is None:
        continue
    a,b,c,cov = fit
    if a >= 0:
        # not a peak (faint) in this window -> unreliable / shape not concave-down
        continue
    t_vertex_offset = -b/(2*a)
    if abs(t_vertex_offset) > HALF_WINDOW:
        continue  # vertex fell outside window -> unreliable extrapolation
    t_obs = t_calc + t_vertex_offset
    oc = t_obs - t_calc

    # propagate uncertainty: t_v = -b/(2a)
    t_obs_err = None
    if cov is not None:
        dtv_da = b/(2*a*a)
        dtv_db = -1.0/(2*a)
        var_tv = (dtv_da**2)*cov['var_a'] + (dtv_db**2)*cov['var_b'] + 2*dtv_da*dtv_db*cov['cov_ab']
        if var_tv > 0:
            t_obs_err = var_tv**0.5

    mag_min = a*t_vertex_offset*t_vertex_offset + b*t_vertex_offset + c

    n_used = len(pts)
    results.append({
        'E': E, 't_calc': t_calc, 't_obs': t_obs, 'OC_days': oc,
        'OC_err_days': t_obs_err, 'n_nights': n_used, 'year': jd_to_decyear(t_obs),
        'mag_min': mag_min
    })

print("Epochs with a fitted primary minimum:", len(results))
for res in results:
    print(res)

# ---- weighted linear fit of O-C vs E: OC = dT0 + dP*E ----
# (a nonzero slope means the trial period needs a small correction;
#  a nonzero intercept means T0 needs a small correction)
we = [1.0/(r['OC_err_days']**2) if r['OC_err_days'] else 1.0 for r in results]
Ee = [r['E'] for r in results]
OCe = [r['OC_days'] for r in results]
Sw = sum(we)
SwE = sum(wi*Ei for wi,Ei in zip(we,Ee))
SwEE = sum(wi*Ei*Ei for wi,Ei in zip(we,Ee))
SwOC = sum(wi*oi for wi,oi in zip(we,OCe))
SwEOC = sum(wi*Ei*oi for wi,Ei,oi in zip(we,Ee,OCe))
denom = Sw*SwEE - SwE*SwE
dP = (Sw*SwEOC - SwE*SwOC)/denom
dT0 = (SwEE*SwOC - SwE*SwEOC)/denom
print(f"\nWeighted linear fit: O-C = {dT0:.4f} + ({dP:.6f})*E days")
print(f"  -> Corrected period estimate: {P+dP:.6f} d")
print(f"  -> Corrected T0 estimate:     {T0+dT0:.4f}")

import csv as csvmod
csv_path = '/Users/mgarrod/CascadeProjects/tcrb/analysis/oc_results.csv'
with open(csv_path, 'w', newline='') as f:
    w = csvmod.writer(f)
    w.writerow(['E','year','t_calc_JD','t_obs_JD','OC_days','OC_err_days','n_nights'])
    for r in results:
        w.writerow([r['E'], f"{r['year']:.4f}", f"{r['t_calc']:.4f}", f"{r['t_obs']:.4f}",
                    f"{r['OC_days']:.4f}", f"{r['OC_err_days']:.4f}" if r['OC_err_days'] else '',
                    r['n_nights']])
print("\nWrote", csv_path)

import json
with open('/Users/mgarrod/CascadeProjects/tcrb/analysis/oc_results.json','w') as f:
    json.dump({'params': {'P': P, 'T0': T0, 'dP': dP, 'dT0': dT0},
                'results': results}, f, indent=2)
