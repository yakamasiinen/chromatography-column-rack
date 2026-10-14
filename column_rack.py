#!/usr/bin/env python3
"""
column_rack.py - parametric wall-hung rack for prepacked chromatography columns
(Cytiva HiTrap-style 1 mL / 5 mL columns), for FDM 3D printing.

    back plate ── hangs on two hooks (e.g. magnetic hooks) through holes near the top
       │
    shelves ── one per row, each with U-shaped slots open to the front; the column hangs
       │       in a slot, and a small lip on the front edge stops it from sliding out
       │
    gussets ── triangular braces under the ends of each shelf
       │
    feet (optional) ── to stand the rack on a bench against a wall instead of hanging it

The rack is built as one sweep along x of a (y, z) cross-section whose points move with x
(slots, lips, gussets, braces) or split around the hook holes, so the mesh is a single closed
surface without any Boolean operations.

Requires Python 3.8+ and numpy (plus matplotlib only for --preview).

Usage
-----
    python column_rack.py                                   # 5 mL, 3 rows x 9  -> column_rack_5mL_3x9.stl
    python column_rack.py --column 5 --rows 2 --cols 7      # 2 rows x 7
    python column_rack.py --column 1                        # 1 mL (slot 10.5 mm, 12 per row)
    python column_rack.py --width 256                       # fit a 256 mm bed (5 mL: 8, 1 mL: 11 per row)
    python column_rack.py --row-types 1,5,5                # mixed: top row 1 mL, two rows 5 mL
    python column_rack.py --set hook_slot=6                 # hook holes as horizontal slots (+6 mm)
    python column_rack.py --feet 2                          # add 2 feet: stands on a bench against a wall
    python column_rack.py --text "Iritani Lab"              # emboss text between the hook holes
    python column_rack.py --list                            # all parameters, auto values resolved
    python column_rack.py --preview rack.png                # also render a preview (matplotlib)

Coordinates (mm): x = along the rack (left to right, seen from the front), y = depth
(y = 0 is the front edge of the shelves, the plate's back face is at y = shelf_depth + plate_t),
z = up (z = 0 is the bottom edge of the plate).
Print lying on its back (back face on the bed): no supports are needed (the lips become short
3 mm overhangs at the top of the print, which print fine without support).
"""
import argparse, struct, sys
import numpy as np

# --------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------
COMMON = dict(
    # ---- layout ----
    column='5',            # column type for every row ('5' or '1'); see COLUMN_TYPES
    row_types=None,        # per-row types, top row first, e.g. ('1','5','5'); overrides column/rows
    rows=3,                # number of shelves
    cols=None,             # slots per shelf. None = as many as fit in `width` (5 mL: 9, 1 mL: 12);
                           # one number: used for the widest-pitch rows, the others are filled to the
                           # same width; or one value per row, top row first
    width=288.0,           # target rack width for the automatic column count (288 = the original rack)
    row_pitch=78.0,        # vertical distance between shelves (or one value per gap, top first)
    # ---- slots, per column type (None = derived / must be given) ----
    slot_w_5=20.0,         # U-slot width for 5 mL (measured from the original rack)
    pitch_5=None,          # slot pitch for 5 mL (auto: slot_w + finger_w = 32)
    slot_w_1=10.5,         # U-slot width for 1 mL: measured OD 10 mm + 0.5 mm clearance
    pitch_1=None,          # slot pitch for 1 mL (auto: slot_w + finger_w)
    finger_w=12.0,         # shelf material between two slots (used for the auto pitch)
    end_extra=0.0,         # extra shelf width beyond the first/last slot, each side
    # ---- shelf ----
    shelf_t=5.0,           # shelf thickness
    shelf_depth=None,      # front edge to plate (auto: back_gap + slot_w + front_run of the widest slot)
    front_run=13.0,        # straight part of the slot from the front edge to the start of the round back
    back_gap=5.0,          # shelf material between the slot back and the plate
    mouth_chamfer=1.5,     # 45 deg lead-in at the slot mouth (0 = sharp corner, as the original)
    lip_h=3.0,             # retaining lip on the front edge of each finger: height (0 = none)
    lip_d=3.0,             # ... and depth
    # ---- strength ----
    brace=3.0,             # 45 deg chamfer under each shelf where it meets the plate (0 = none)
    top_chamfer=1.0,       # 45 deg chamfer on top of each shelf at the plate (0 = none)
    gussets='ends',        # triangular braces under the shelves: 'ends', 'all' (every finger) or 'none'
    gusset_h=30.0,         # gusset height at the plate
    # ---- back plate ----
    plate_t=5.0,           # plate thickness
    plate_below=None,      # plate below the lowest shelf (auto: gusset_h, or 10 without gussets;
                           # with feet: room for the bottom row of columns, see foot_clear)
    # ---- feet (optional): stand the rack on a bench against a wall ----
    feet=0,                # number of feet (0 = wall-hung only); spread evenly, the outer ones at the ends
    foot_w=12.0,           # foot width (x); inner feet sit under a finger between two slots
    foot_t=5.0,            # foot thickness (z)
    foot_ext=20.0,         # how far the feet reach in front of the shelf front edge
    foot_brace=15.0,       # 45 deg brace between foot and plate (size along the plate)
    foot_clear=None,       # with feet: free height under the lowest shelf (auto: col_below + 5, else 50)
    plate_above=20.0,      # plate above the top of the top shelf
    # ---- hook holes ----
    hook_n=2,              # number of holes (0 = none)
    hook_d=8.0,            # hole diameter (height)
    hook_slot=0.0,         # extra horizontal length: > 0 makes slots, which tolerate inaccurate hook spacing
    hook_from_top=10.0,    # hole centre below the top edge of the plate
    hook_inset=None,       # hole centre from the plate end (auto: above the first inner finger, as the original)
    # ---- embossed text (optional, needs matplotlib) ----
    text='',               # e.g. 'Iritani Lab'; '' = none. Goes on the plate front between the hook holes
    text_size=None,        # font size in mm (auto: up to 10 mm, smaller if needed to fit between the holes)
    text_h=0.6,            # emboss height
    text_z=None,           # height of the text centre (auto: level with the hook holes)
    text_font='DejaVu Sans:bold',  # 'family:weight'; for Japanese e.g. 'Yu Gothic:bold' or 'Meiryo:bold'
    # ---- checks (optional) ----
    bed_x=256.0, bed_y=256.0,   # printer bed (the rack is printed lying on its back: W x H)
    col_mass=None,         # mass of one filled column (g) for the load estimate
    col_above=None,        # how far a hanging column sticks up above the shelf top (mm)
    col_below=None,        # how far it hangs below the shelf bottom (mm)
    shoulder_d=None,       # diameter of the part that rests on the shelf (must clear neighbours/plate)
    density=1.24,          # filament density g/cm3 (PLA 1.24, PETG 1.27)
    # ---- mesh ----
    arc_seg=24,            # segments per slot back (half circle)
    hole_seg=24,           # segments per hook hole (full circle)
)

# Column types. 5 mL values were measured from the original SketchUp rack (column_rack.stl):
# U-slot 20 mm wide, 32 mm pitch, centre 23 mm behind the front edge.
COLUMN_TYPES = ('5', '1')
AUTO_KEYS = ('cols', 'pitch_5', 'pitch_1', 'shelf_depth', 'plate_below', 'hook_inset', 'foot_clear')
EPS = 0.005            # x offset used where the cross-section changes abruptly (slot walls)


def as_list(v, n, name):
    if isinstance(v, (tuple, list)):
        if len(v) != n: sys.exit(f'{name} needs {n} values (one per row/gap), got {len(v)}')
        return [float(x) for x in v]
    return [float(v)] * n


def resolve(p):
    """Fill in derived values and the per-row layout."""
    types = [str(t).replace('.0', '') for t in (p['row_types'] or [p['column']] * int(p['rows']))]
    for t in types:
        if t not in COLUMN_TYPES: sys.exit(f'unknown column type {t!r} (use one of {COLUMN_TYPES})')
        if p[f'slot_w_{t}'] is None:
            sys.exit(f'slot width for {t} mL columns is not set. Measure the part of the column that sits '
                     f'in the slot and give it with --slot-w (or --set slot_w_{t}=VALUE); see README.')
        if p[f'pitch_{t}'] is None:
            p[f'pitch_{t}'] = round(p[f'slot_w_{t}'] + p['finger_w'], 2)
    p['rows'] = len(types)
    if p['cols'] is None:   # as many slots as fit in the target width
        cols = [max(1, int(np.floor((p['width'] - 2 * p['end_extra']) / p[f'pitch_{t}'] + 1e-9))) for t in types]
        p['cols'] = cols[0] if len(set(cols)) == 1 else tuple(cols)
    elif isinstance(p['cols'], (tuple, list)):
        cols = [int(c) for c in as_list(p['cols'], len(types), 'cols')]
    else:   # one number: it applies to the widest-pitch type; other rows are filled to the same width
        wmax = int(p['cols']) * max(p[f'pitch_{t}'] for t in types)
        cols = [int(p['cols']) if p[f'pitch_{t}'] == max(p[f'pitch_{u}'] for u in types)
                else int(np.floor(wmax / p[f'pitch_{t}'] + 1e-9)) for t in types]
    gaps = as_list(p['row_pitch'], max(len(types) - 1, 1), 'row_pitch')[:len(types) - 1]
    if p['shelf_depth'] is None:
        p['shelf_depth'] = round(p['back_gap'] + max(p[f'slot_w_{t}'] for t in types) + p['front_run'], 2)
    below = p['gusset_h'] if p['gussets'] != 'none' else p['brace']
    if p['feet'] and p['foot_clear'] is None:
        p['foot_clear'] = round(p['col_below'] + 5, 1) if p['col_below'] else 50.0
    if p['plate_below'] is None:
        if p['feet']:   # the lowest gusset/brace must end above the foot brace
            p['plate_below'] = round(max(p['foot_clear'] + p['foot_t'],
                                         below + p['foot_t'] + p['foot_brace'] + 2), 1)
        else:
            p['plate_below'] = p['gusset_h'] if p['gussets'] != 'none' else 10.0
    rows = []
    W = max(c * p[f'pitch_{t}'] for t, c in zip(types, cols)) + 2 * p['end_extra']
    # rows are listed top first; z is built bottom up
    zb = p['plate_below']
    zs = []
    for i in range(len(types)):
        zs.append(zb)
        if i < len(gaps): zb += gaps[::-1][i]
    zs = zs[::-1]
    for t, c, z in zip(types, cols, zs):
        pitch, w = p[f'pitch_{t}'], p[f'slot_w_{t}']
        x0 = (W - c * pitch) / 2
        r = w / 2
        rows.append(dict(type=t, n=c, pitch=pitch, w=w, r=r, x0=x0,
                         centres=[x0 + pitch * (i + 0.5) for i in range(c)],
                         yc=p['shelf_depth'] - p['back_gap'] - r,       # centre of the round slot back
                         zb=z, zt=z + p['shelf_t']))
    p['_rows'] = rows
    p['_W'] = W
    p['_H'] = rows[0]['zt'] + p['plate_above']
    if p['hook_inset'] is None:
        top = rows[0]
        p['hook_inset'] = round(top['x0'] + top['pitch'], 2) if top['n'] > 2 else round(W / 4, 2)
    return p


def validate(p):
    err = []
    for k in ('shelf_t', 'plate_t', 'back_gap', 'front_run', 'finger_w'):
        if p[k] <= 0: err.append(f'{k} must be > 0')
    for k in ('mouth_chamfer', 'lip_h', 'lip_d', 'brace', 'top_chamfer', 'gusset_h', 'hook_slot', 'end_extra'):
        if p[k] < 0: err.append(f'{k} must be >= 0')
    if p['gussets'] not in ('ends', 'all', 'none'): err.append("gussets must be 'ends', 'all' or 'none'")
    if err: sys.exit('invalid parameters: ' + '; '.join(err))


def check_layout(p):
    """Geometric problems that would make an invalid cross-section (list of strings)."""
    probs = []; rows = p['_rows']; YF = p['shelf_depth']
    for i, R in enumerate(rows):
        tag = f'row {i + 1} ({R["type"]} mL)'
        if R['pitch'] - R['w'] - 2 * p['mouth_chamfer'] < 1.0:
            probs.append(f'{tag}: fingers between slots are too narrow (pitch - slot_w - 2*mouth_chamfer < 1 mm)')
        if R['yc'] < p['mouth_chamfer'] + 0.5:
            probs.append(f'{tag}: slot is deeper than the shelf (increase shelf_depth or front_run)')
        back = R['yc'] + R['r']
        if back + max(p['brace'], p['top_chamfer']) + 0.3 > YF:
            probs.append(f'{tag}: back_gap ({p["back_gap"]}) must be larger than brace/top_chamfer')
        if p['lip_d'] + p['mouth_chamfer'] >= min(YF - p['top_chamfer'] - 1, R['yc'] - 0.5):
            probs.append(f'{tag}: lip_d (+ mouth_chamfer) reaches into the round back of the slot')
        if p['gussets'] != 'none' and p['gusset_h'] < p['brace']:
            probs.append('gusset_h must be >= brace')
    for i in range(len(rows) - 1):
        up, lo = rows[i], rows[i + 1]
        below = p['gusset_h'] if p['gussets'] != 'none' else p['brace']
        if up['zb'] - below <= lo['zt'] + max(p['lip_h'], p['top_chamfer']) + 1:
            probs.append(f'rows {i + 1} and {i + 2} are too close (row_pitch) for the gussets/lips')
    need = p['gusset_h'] if p['gussets'] != 'none' else p['brace']
    if p['feet']:
        need += p['foot_t'] + p['foot_brace'] + 1
        if p['foot_ext'] + p['shelf_depth'] < p['foot_brace'] + 1: probs.append('feet too short (foot_ext)')
        for a, b in foot_spans(p):
            if b - a < 2: probs.append('foot_w too small')
        fs = foot_spans(p)
        if any(b2[0] < a2[1] for a2, b2 in zip(fs, fs[1:])): probs.append('feet overlap: fewer feet or smaller foot_w')
    if p['plate_below'] < need - 1e-9: probs.append(f'plate_below must be >= {need:g} (gusset_h/brace, + foot_t + foot_brace with feet)')
    if p['hook_n']:
        zh = p['_H'] - p['hook_from_top']
        if zh + p['hook_d'] / 2 > p['_H'] - 2: probs.append('hook holes too close to the top edge (hook_from_top)')
        if zh - p['hook_d'] / 2 < rows[0]['zt'] + p['top_chamfer'] + 2:
            probs.append('hook holes overlap the top shelf: increase plate_above or reduce hook_from_top')
        for x0, x1 in hole_spans(p):
            if x0 < 2 or x1 > p['_W'] - 2: probs.append('hook holes too close to the plate ends (hook_inset)')
        sp = hole_spans(p)
        for a, b in zip(sp, sp[1:]):
            if b[0] - a[1] < 3: probs.append('hook holes overlap each other')
    return probs

# --------------------------------------------------------------------------------------
# Cross-section along x
# --------------------------------------------------------------------------------------
def hole_outline(p):
    """Vertices (x, z) of the upper half of a hook hole centred at x = 0: a polygonal circle with
    vertical flats at its left/right ends (so the cross-section stays valid), stretched by hook_slot."""
    n = max(8, int(p['hole_seg']) // 2 * 2)
    R = p['hook_d'] / 2 / np.cos(np.pi / n)            # circumradius: the flats are at hook_d/2
    th = (np.arange(n // 2) + 0.5) * 2 * np.pi / n      # upper half, from right to left
    x = R * np.cos(th) + np.where(np.cos(th) > 0, 1, -1) * p['hook_slot'] / 2
    z = R * np.sin(th)
    return x[::-1], z[::-1]                            # left to right


def hole_centres(p):
    if not p['hook_n']: return []
    n = int(p['hook_n']); W = p['_W']; a = p['hook_inset']
    if n == 1: return [W / 2]
    return list(np.linspace(a, W - a, n))


def hole_spans(p):
    hx, _ = hole_outline(p)
    return [(c + hx[0], c + hx[-1]) for c in hole_centres(p)]


def foot_spans(p):
    """x ranges of the feet (outer feet flush with the plate ends)."""
    n = int(p['feet']); W = p['_W']; w = min(p['foot_w'], W)
    if n <= 0: return []
    cs = [W / 2] if n == 1 else list(np.linspace(w / 2, W - w / 2, n))
    # inner feet go under a finger of the lowest row (between two slots), not under a hanging column
    R = p['_rows'][-1]
    fingers = [R['x0'] + R['pitch'] * k for k in range(1, R['n'])]
    if fingers:
        for i in range(len(cs)):
            if (n == 1 or 0 < i < n - 1) and w < p['_W']:
                cs[i] = min(fingers, key=lambda f: abs(f - cs[i]))
    return [(c - w / 2, c + w / 2) for c in cs]


def foot_on(p, x):
    W = p['_W']
    for a, b in foot_spans(p):
        a = -1 if a <= 1e-9 else a; b = W + 1 if b >= W - 1e-9 else b
        if a + EPS / 2 <= x <= b - EPS / 2: return True
    return False


def stations(p):
    """x positions of all cross-sections."""
    W = p['_W']; X = [0.0, W]; mc = p['mouth_chamfer']
    for R in p['_rows']:
        n_arc = max(4, int(p['arc_seg']))
        for c in R['centres']:
            h = R['r']
            X += [c - h - mc, c + h + mc, c - h, c - h + EPS / 2, c - h + EPS, c + h - EPS, c + h - EPS / 2, c + h]
            th = np.linspace(0, np.pi, n_arc + 1)[1:-1]
            X += list(c - h * np.cos(th))
    hx, _ = hole_outline(p)
    for c in hole_centres(p): X += list(c + hx)
    for a, b in foot_spans(p):            # off at a, on at a + EPS (a foot appears in one thin step)
        X += [a, a + EPS, b - EPS, b]
    X = np.unique(np.round(np.clip(X, 0, W), 6))
    return X


def row_front(R, x, p):
    """(y of the shelf front at x, state, on an end finger) for one row.
    state: 'finger', 'slot', or 'mid' = the thin step at a slot wall where the lip has already
    collapsed but the front has not yet moved back (keeps the surface from folding over itself)."""
    mc = p['mouth_chamfer']; h = R['r']
    i = int(np.argmin([abs(x - c) for c in R['centres']]))
    d = abs(x - R['centres'][i])
    end = x <= R['centres'][0] or x >= R['centres'][-1]
    if d > h - 0.25 * EPS:                              # finger (with the 45 deg mouth chamfer)
        return max(0.0, h + mc - d), 'finger', end
    if d > h - 0.75 * EPS:
        return mc, 'mid', end
    if d > h - 1.5 * EPS:                               # first point of the slot: straight wall
        return R['yc'], 'slot', False
    return R['yc'] + np.sqrt(max(h * h - d * d, 0.0)), 'slot', False


def hole_band(p, x):
    """(z1, z2, inside) of the hook-hole band at x."""
    hx, hz = hole_outline(p)
    zh = p['_H'] - p['hook_from_top']
    for c in hole_centres(p):
        if c + hx[0] - 1e-9 <= x <= c + hx[-1] + 1e-9:
            dz = np.interp(x - c, hx, hz)
            return zh - dz, zh + dz, True
    return zh - hz[0], zh + hz[0], False


def section(p, x):
    """Cross-section at x: array of (y, z) points, in the same order at every x (CCW in y-z).
    Point order: back bottom, back hole points, back top, front top, front hole points,
    then for each row (top first) 8 points, then front bottom."""
    YF = p['shelf_depth']; YB = YF + p['plate_t']; H = p['_H']
    z1, z2, _ = hole_band(p, x) if p['hook_n'] else (H - 1, H - 0.5, False)
    pts = [(YB, 0.0), (YB, z1), (YB, z2), (YB, H), (YF, H), (YF, z2), (YF, z1)]
    tc, br, ld, lh = p['top_chamfer'], p['brace'], p['lip_d'], p['lip_h']
    for R in p['_rows']:
        yf, state, end = row_front(R, x, p)
        zt, zb = R['zt'], R['zb']
        gus = state != 'slot' and (p['gussets'] == 'all' or (p['gussets'] == 'ends' and end))
        pts += [(YF, zt + tc), (YF - tc, zt)]                       # top chamfer at the plate
        if state == 'slot' or lh == 0:
            pts += [(yf, zt)] * 3                                    # no lip (collapsed points)
        elif state == 'mid':
            pts += [(yf + ld, zt), (yf + ld, zt), (yf, zt)]          # lip flattened onto the shelf top
        else:
            pts += [(yf + ld, zt), (yf + ld, zt + lh), (yf, zt + lh)]
        pts += [(yf, zb)]                                            # front bottom edge
        if gus and p['gusset_h'] > 0:
            pts += [(yf, zb), (YF, zb - p['gusset_h'])]              # gusset: straight down to the plate
        else:
            pts += [(YF - br, zb), (YF, zb - br)]                    # 45 deg brace under the shelf
    # foot: brace down the plate, top of the foot, its front end; collapsed to (YF, 0) where there is none
    if foot_on(p, x):
        ft, fb, yff = p['foot_t'], p['foot_brace'], -p['foot_ext']
        pts += [(YF, ft + fb), (YF - fb, ft), (yff, ft), (yff, 0.0)]
    else:
        pts += [(YF, 0.0)] * 4
    pts += [(YF, 0.0)]
    return np.array(pts, float)

# --------------------------------------------------------------------------------------
# Mesh
# --------------------------------------------------------------------------------------
def ear_clip(P):
    P = np.asarray(P, float); idx = list(range(len(P))); T = []
    cr = lambda o, a, b: (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    eps = 1e-10
    def blocks(a, b, c, k):
        A, B, C, K = P[a], P[b], P[c], P[k]
        if np.allclose(K, A) or np.allclose(K, C): return False
        c1, c2, c3 = cr(A, B, K), cr(B, C, K), cr(C, A, K)
        if c1 > eps and c2 > eps and c3 > eps: return True
        if abs(c3) <= eps and min(A[0], C[0]) - eps <= K[0] <= max(A[0], C[0]) + eps \
                and min(A[1], C[1]) - eps <= K[1] <= max(A[1], C[1]) + eps: return True
        return False
    while len(idx) > 3:
        for j in range(len(idx)):
            a, b, c = idx[j-1], idx[j], idx[(j+1) % len(idx)]
            if cr(P[a], P[b], P[c]) <= eps: continue
            if any(blocks(a, b, c, k) for k in idx if k not in (a, b, c)): continue
            T.append((a, b, c)); idx.pop(j); break
        else:
            raise RuntimeError('cross-section is not a simple polygon (check parameters)')
    T.append(tuple(idx)); return T


def cap(P2, x, flip):
    """Triangulated end cap of a cross-section at x (duplicate points removed first)."""
    keep = [i for i in range(len(P2)) if np.hypot(*(P2[i] - P2[i - 1])) > 1e-9]
    Q = P2[keep]
    tris = np.array(ear_clip(Q))
    T = np.concatenate([np.full((len(Q), 1), x), Q], 1)[tris]
    return T[:, ::-1] if flip else T


def build_mesh(p):
    X = stations(p)
    S = np.array([section(p, x) for x in X])                         # (n_x, m, 2)
    m = S.shape[1]
    V = np.concatenate([np.repeat(X[:, None, None], m, 1), S], 2)   # (n_x, m, 3)
    eA = [(i, (i + 1) % m) for i in range(m)]
    # inside a hole the back/front hole edges are replaced by the hole's bottom and top surfaces
    eB = [e for e in eA if e not in ((1, 2), (5, 6))] + [(1, 6), (5, 2)]
    spans = hole_spans(p)
    tris = []
    for k in range(len(X) - 1):
        xm = (X[k] + X[k + 1]) / 2
        E = eB if any(a < xm < b for a, b in spans) else eA
        for i, j in E:
            a, b, c, d = V[k, i], V[k, j], V[k + 1, j], V[k + 1, i]
            tris += [(a, b, c), (a, c, d)]
    tris = np.array(tris)
    caps = [cap(S[0], X[0], True), cap(S[-1], X[-1], False)]
    for a, b in spans:                                               # vertical flats at the hole ends
        for xx, flip in ((a, False), (b, True)):
            k = int(np.argmin(np.abs(X - xx))); q = V[k]
            quad = np.array([(q[1], q[2], q[5]), (q[1], q[5], q[6])])
            caps.append(quad[:, ::-1] if flip else quad)
    T = np.concatenate([tris] + caps)
    # weld identical points, drop collapsed triangles
    flat = np.round(T.reshape(-1, 3), 6)
    U, inv = np.unique(flat, axis=0, return_inverse=True)
    F = inv.reshape(-1, 3)
    ok = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])
    F = F[ok]
    T = U[F]
    vol = np.einsum('ij,ij->i', T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6
    if vol < 0: F = F[:, ::-1]; T = U[F]; vol = -vol
    return U, F, T, vol


def check_manifold(F):
    E = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]).astype(np.int64)
    n = F.max() + 1
    key = E[:, 0] * n + E[:, 1]; rkey = E[:, 1] * n + E[:, 0]
    u, cnt = np.unique(key, return_counts=True)
    return bool((cnt == 1).all() and np.isin(rkey, u).all())


def sections_simple(p):
    """True if no cross-section along x intersects itself (zero-length edges are ignored)."""
    for x in stations(p):
        P = section(p, x)
        keep = [i for i in range(len(P)) if np.hypot(*(P[i] - P[i - 1])) > 1e-9]
        Q = P[keep]; A = Q; B = np.roll(Q, -1, 0); n = len(Q)
        cr = lambda o, a, b: (a[..., 0]-o[..., 0])*(b[..., 1]-o[..., 1]) - (a[..., 1]-o[..., 1])*(b[..., 0]-o[..., 0])
        a, b = A[:, None], B[:, None]; c, d = A[None, :], B[None, :]
        hit = (cr(a, b, c) * cr(a, b, d) < -1e-12) & (cr(c, d, a) * cr(c, d, b) < -1e-12)
        i, j = np.nonzero(np.triu(hit, 2))
        if len(i) and not all((jj + 1) % n == ii for ii, jj in zip(i, j)): return False
    return True


def sliver_count(T, tol=1e-6):
    A = np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1) / 2
    return int((A < tol).sum())

# --------------------------------------------------------------------------------------
# Optional embossed text (same method as tube_stand.py)
# --------------------------------------------------------------------------------------
def poly_area(P):
    P = np.asarray(P, float)
    return 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])


def _pip(pt, poly):
    x, y = pt; inside = False; n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]; x1, y1 = poly[(i + 1) % n]
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0): inside = not inside
    return inside

def _seg_cross(a, b, c, d):
    cr = lambda o, p, q: (p[0]-o[0])*(q[1]-o[1]) - (p[1]-o[1])*(q[0]-o[0])
    return cr(a, b, c) * cr(a, b, d) < 0 and cr(c, d, a) * cr(c, d, b) < 0

def _merge_holes(outer, holes):
    """Bridge the holes of a glyph (e.g. inside 'o') into its outline so it can be ear-clipped."""
    poly = [tuple(q) for q in outer]
    for hole in sorted(holes, key=lambda h: -max(q[0] for q in h)):
        hole = [tuple(q) for q in hole]; mi = max(range(len(hole)), key=lambda i: hole[i][0]); M = hole[mi]
        edges = [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]
        for h in holes: edges += [(tuple(h[i]), tuple(h[(i + 1) % len(h)])) for i in range(len(h))]
        counts = {}
        for q in poly: counts[q] = counts.get(q, 0) + 1
        for vi, V in sorted(enumerate(poly), key=lambda t: (t[1][0] - M[0])**2 + (t[1][1] - M[1])**2):
            if counts[V] > 1: continue                 # never reuse an earlier bridge point
            if all(not _seg_cross(M, V, e0, e1) for e0, e1 in edges if V not in (e0, e1) and M not in (e0, e1)): break
        else:
            raise RuntimeError('could not bridge a glyph hole')
        poly = poly[:vi + 1] + hole[mi:] + hole[:mi + 1] + [V] + poly[vi + 1:]
    return poly

def glyph_outlines(text, size, font):
    """Text outlines (mm) from matplotlib's TextPath: list of (outer, [holes])."""
    try:
        from matplotlib.textpath import TextPath
        from matplotlib.font_manager import FontProperties
        from matplotlib.path import Path as MPath
    except ImportError:
        sys.exit('--text needs matplotlib (pip install matplotlib)')
    fam, _, wt = font.partition(':')
    prop = FontProperties(family=fam, weight=wt or 'normal')
    # refuse silently-substituted fonts and missing characters (they would print as boxes)
    import matplotlib.font_manager as fm
    from matplotlib.ft2font import FT2Font
    try:
        path = fm.findfont(prop, fallback_to_default=False)
    except Exception:
        sys.exit(f'font "{fam}" not found. Installed families include: '
                 + ', '.join(sorted({f.name for f in fm.fontManager.ttflist})[:15]) + ' ...')
    ft = FT2Font(path)
    missing = sorted({ch for ch in text if not ch.isspace() and ft.get_char_index(ord(ch)) == 0})
    if missing:
        sys.exit(f'font "{fam}" has no glyphs for: {"".join(missing)}. Choose a font that has them, '
                 f'for Japanese e.g. --set "text_font=Yu Gothic:bold" (Windows), "Hiragino Sans" (macOS) '
                 f'or "Noto Sans CJK JP:bold" (Linux).')
    tp = TextPath((0, 0), text, size=size, prop=prop)
    polys, cur = [], []
    for seg, code in tp.iter_bezier():                 # flatten curves finely
        if code == MPath.MOVETO:
            if len(cur) > 3: polys.append(np.array(cur))
            cur = [tuple(seg.control_points[0])]; continue
        if code == MPath.CLOSEPOLY: continue
        ts = np.linspace(0, 1, 2 if seg.degree == 1 else 9)[1:]
        cur += [tuple(q) for q in np.atleast_2d(seg(ts))]
    if len(cur) > 3: polys.append(np.array(cur))
    polys = [q[:-1] if np.allclose(q[0], q[-1]) else q for q in polys]
    polys = [q[np.r_[True, np.any(np.abs(np.diff(q, axis=0)) > 1e-6, axis=1)]] for q in polys]
    # drop degenerate contours (some fonts, e.g. Yu Gothic, contain tiny or zero-area pieces)
    polys = [q for q in polys if len(q) >= 3 and abs(poly_area(q)) > 1e-4]
    if not polys: sys.exit('text produced no outlines')
    # filled contours vs holes by winding direction (fonts use the non-zero rule: holes run the
    # other way). Overlapping filled contours (common in CJK fonts) become separate overlapping shells.
    # Decide per contour, from the largest down: a contour inside nothing is filled; a contour
    # inside another one has the same role as its container if it winds the same way, else the
    # opposite role. Works for both winding conventions, for fallback glyphs from other fonts,
    # and for overlapping strokes (CJK), which become separate overlapping shells.
    areas = [poly_area(q) for q in polys]
    order = sorted(range(len(polys)), key=lambda i: -abs(areas[i]))
    inside = lambda i, j: all(_pip(q, polys[j]) for q in polys[i][::max(1, len(polys[i]) // 8)])
    filled_role, container = {}, {}
    for k, i in enumerate(order):
        cands = [j for j in order[:k] if abs(areas[j]) > abs(areas[i]) and inside(i, j)]
        if not cands: filled_role[i] = True; continue
        c = min(cands, key=lambda j: abs(areas[j])); container[i] = c
        filled_role[i] = filled_role[c] if np.sign(areas[i]) == np.sign(areas[c]) else not filled_role[c]
    polys = [q if filled_role[i] == (areas[i] > 0) else q[::-1] for i, q in enumerate(polys)]  # filled CCW, holes CW
    groups = {i: [] for i in range(len(polys)) if filled_role[i]}
    for i in range(len(polys)):
        if not filled_role[i] and i in container and filled_role[container[i]]: groups[container[i]].append(i)
    allp = np.concatenate(polys)
    return [(polys[o], [polys[h] for h in hs]) for o, hs in groups.items()], allp

def text_area(p):
    """(x0, x1, z0, z1, zc) of the free band on the plate front for the text."""
    W, H = p['_W'], p['_H']
    R = p['_rows'][0]
    z0 = R['zt'] + p['top_chamfer'] + 1.0
    z1 = H - 1.0
    sp = hole_spans(p)
    if len(sp) >= 2:   # between the two middle-most holes around the centre
        left = max([b for a, b in sp if b <= W / 2 + 1e-9], default=2.0)
        right = min([a for a, b in sp if a >= W / 2 - 1e-9], default=W - 2.0)
        x0, x1 = left + 3.0, right - 3.0
    elif len(sp) == 1:
        x0, x1 = 3.0, sp[0][0] - 3.0
    else:
        x0, x1 = 3.0, W - 3.0
    zc = p['text_z'] if p['text_z'] is not None else \
        (H - p['hook_from_top'] if p['hook_n'] else (z0 + z1) / 2)
    return x0, x1, z0, z1, zc


def text_triangles(p):
    """Embossed text on the plate front, as closed shells sunk 0.4 mm into the plate
    (the slicer merges them). Returns (triangles, info)."""
    x0, x1, z0, z1, zc = text_area(p)
    size = p['text_size']
    auto = size is None
    if auto:
        size = 10.0
        g, allp = glyph_outlines(p['text'], size, p['text_font'])
        wid = allp[:, 0].max() - allp[:, 0].min(); tall = allp[:, 1].max() - allp[:, 1].min()
        f = min(1.0, (x1 - x0) / wid, (2 * min(zc - z0, z1 - zc)) / tall)
        if f < 1: size = max(3.0, np.floor(size * f * 10) / 10)
    glyphs, allp = glyph_outlines(p['text'], size, p['text_font'])
    uc = (allp[:, 0].min() + allp[:, 0].max()) / 2; vc = (allp[:, 1].min() + allp[:, 1].max()) / 2
    xc = (x0 + x1) / 2 if len(hole_spans(p)) != 1 else (x0 + x1) / 2
    YF = p['shelf_depth']; sink = 0.4
    def to3d(uv, off):   # seen from the front: u -> +x, v -> +z, emboss towards -y
        return np.stack([xc + uv[:, 0] - uc, np.full(len(uv), YF - off), zc + uv[:, 1] - vc], -1)
    tris = []
    for outer, holes in glyphs:
        merged = np.array(_merge_holes(outer, holes))
        if len(merged) < 3: continue
        capi = np.array([t for t in ear_clip(merged) if len(t) == 3])
        if len(capi) == 0: continue
        top, bot = to3d(merged, p['text_h']), to3d(merged, -sink)
        tris += [top[capi], bot[capi][:, ::-1]]
        for ring in [outer] + holes:
            b_, t_ = to3d(ring, -sink), to3d(ring, p['text_h'])
            b2, t2 = np.roll(b_, -1, 0), np.roll(t_, -1, 0)
            tris += [np.stack([b_, b2, t2], 1), np.stack([b_, t2, t_], 1)]
    T = np.concatenate(tris)
    if np.einsum('ij,ij->i', T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() < 0: T = T[:, ::-1]
    info = dict(size=size, auto=auto, x=(T[..., 0].min(), T[..., 0].max()), z=(T[..., 2].min(), T[..., 2].max()),
                area=(x0, x1, z0, z1))
    return T, info


def check_text(info):
    probs = []
    x0, x1, z0, z1 = info['area']
    if info['x'][0] < x0 - 1e-6 or info['x'][1] > x1 + 1e-6:
        probs.append(f'text is {info["x"][1] - info["x"][0]:.0f} mm wide but only {x1 - x0:.0f} mm is free between '
                     f'the hook holes: shorten it, reduce text_size, or move the holes (hook_inset)')
    if info['z'][0] < z0 - 1e-6 or info['z'][1] > z1 + 1e-6:
        probs.append(f'text (z {info["z"][0]:.1f}-{info["z"][1]:.1f}) does not fit between the top shelf and the top '
                     f'edge (z {z0:.1f}-{z1:.1f}): reduce text_size, change text_z or increase plate_above')
    if info['size'] < 4:
        probs.append(f'text size {info["size"]:g} mm is small for a 0.4 mm nozzle (4 mm or more is easier to read)')
    return probs

# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------
def write_stl(path, T, name='column rack'):
    T = np.asarray(T, np.float32)
    nrm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); L = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = nrm / np.where(L == 0, 1, L)
    rec = np.zeros(len(T), np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]))
    rec['n'] = nrm; rec['v'] = T
    with open(path, 'wb') as f:
        f.write(name.encode()[:80].ljust(80, b' ')); f.write(struct.pack('<I', len(T))); f.write(rec.tobytes())


def render(T, az, el, width=900, light=(-0.45, 0.55, 0.7)):
    """Tiny orthographic z-buffer renderer (numpy only). Returns an RGB float image."""
    az, el = np.radians(az), np.radians(el)
    # camera: looks along -v; x' to the right, y' up
    v = np.array([np.sin(az) * np.cos(el), -np.cos(az) * np.cos(el), np.sin(el)])
    r = np.cross([0, 0, 1], v); r /= np.linalg.norm(r); u = np.cross(v, r)
    P = T.reshape(-1, 3) - T.reshape(-1, 3).mean(0)
    sx, sy, sz = P @ r, P @ u, P @ v
    s = (width - 20) / (sx.max() - sx.min()); h = int((sy.max() - sy.min()) * s) + 20
    X = ((sx - sx.min()) * s + 10).reshape(-1, 3); Y = ((sy.max() - sy) * s + 10).reshape(-1, 3)
    Z = sz.reshape(-1, 3)
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    Lw = np.array(light, float); Lw /= np.linalg.norm(Lw)
    L = Lw[0] * r + Lw[1] * u + Lw[2] * v                      # light fixed relative to the camera
    shade = 0.28 + 0.72 * np.abs(n @ L)
    zbuf = np.full((h, width), -np.inf); img = np.ones((h, width, 3))
    base = np.array([0.93, 0.86, 0.70])
    for k in np.nonzero(np.abs(n @ v) > 1e-9)[0]:
        x, y, z = X[k], Y[k], Z[k]
        x0, x1 = int(max(np.floor(x.min()), 0)), int(min(np.ceil(x.max()), width - 1))
        y0, y1 = int(max(np.floor(y.min()), 0)), int(min(np.ceil(y.max()), h - 1))
        if x1 < x0 or y1 < y0: continue
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
        d = (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2])
        if abs(d) < 1e-12: continue
        w0 = ((y[1] - y[2]) * (gx - x[2]) + (x[2] - x[1]) * (gy - y[2])) / d
        w1 = ((y[2] - y[0]) * (gx - x[2]) + (x[0] - x[2]) * (gy - y[2])) / d
        w2 = 1 - w0 - w1
        m = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
        if not m.any(): continue
        zz = w0 * z[0] + w1 * z[1] + w2 * z[2]
        sub = zbuf[y0:y1 + 1, x0:x1 + 1]
        m &= zz > sub
        sub[m] = zz[m]
        img[y0:y1 + 1, x0:x1 + 1][m] = base * shade[k]
    return img


def preview(path, T, p):
    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit('--preview needs matplotlib (pip install matplotlib)')
    # model: front edge at y = 0, plate behind it (+y). View from the front, slightly from above.
    views = [(-28, 22, 'front left'), (24, 30, 'front right, from above')]
    imgs = [render(T, az, el) for az, el, _ in views]
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.2))
    for ax, im, (_, _, t) in zip(axs, imgs, views):
        ax.imshow(im); ax.set_axis_off(); ax.set_title(t, fontsize=10)
    plt.tight_layout(); plt.savefig(path, dpi=110); plt.close(fig)

# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------
def parse_value(s):
    if s.strip().lower() in ('none', 'null', ''): return None
    if ',' in s:
        parts = [x.strip() for x in s.split(',')]
        try: return tuple(float(x) for x in parts)
        except ValueError: return tuple(parts)
    try: return float(s)
    except ValueError: return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--column', choices=COLUMN_TYPES, help='column type for all rows (mL)')
    ap.add_argument('--rows', type=int, help='number of shelves')
    ap.add_argument('--feet', type=int, help='add N feet to stand the rack against a wall (0 = wall-hung only)')
    ap.add_argument('--cols', type=str, help='slots per shelf (or comma list per row, top first); '
                                             'default: as many as fit in --width')
    ap.add_argument('--width', type=float, help='target width for the automatic column count (default 288 mm)')
    ap.add_argument('--row-types', type=str, help='column type per row, top first, e.g. 1,5,5')
    ap.add_argument('--slot-w', type=float, help='slot width (mm) for the column type(s) used')
    ap.add_argument('--text', default=None, help='emboss this text between the hook holes (needs matplotlib)')
    ap.add_argument('--set', action='append', default=[], metavar='KEY=VALUE', help='override any parameter')
    ap.add_argument('-o', '--out', default=None, help='output STL path')
    ap.add_argument('--preview', default=None, help='also write a PNG preview (needs matplotlib)')
    ap.add_argument('--profile-csv', default=None, help='write the cross-section at a finger and a slot centre')
    ap.add_argument('--list', action='store_true', help='print all parameters (auto values resolved) and exit')
    a = ap.parse_args(argv)
    p = dict(COMMON)
    if a.column: p['column'] = a.column
    if a.rows is not None: p['rows'] = a.rows
    if a.feet is not None: p['feet'] = a.feet
    if a.cols is not None: p['cols'] = parse_value(a.cols)
    if a.width is not None: p['width'] = a.width
    if a.text is not None: p['text'] = a.text
    if a.row_types: p['row_types'] = tuple(x.strip() for x in a.row_types.split(','))
    for kv in a.set:
        if '=' not in kv: sys.exit(f'--set expects KEY=VALUE, got: {kv}')
        k_, v_ = kv.split('=', 1)
        if k_ not in p: sys.exit(f'unknown parameter: {k_}')
        p[k_] = v_.strip() if isinstance(COMMON[k_], str) else parse_value(v_)   # text parameters stay strings
    if isinstance(p['row_types'], (str, float)): p['row_types'] = tuple(str(p['row_types']).split(','))
    if p['row_types']: p['row_types'] = tuple(str(t).replace('.0', '') for t in p['row_types'])
    p['column'] = str(p['column']).replace('.0', '')
    if a.slot_w is not None:
        for t in set(p['row_types'] or [p['column']]): p[f'slot_w_{t}'] = a.slot_w
    for k_ in ('hook_n', 'arc_seg', 'hole_seg', 'rows', 'feet'): p[k_] = int(p[k_])
    validate(p)
    resolve(p)
    if a.list:
        fmt = lambda v: f'{v:g}' if isinstance(v, float) else v
        for k_ in COMMON:
            auto = k_ in AUTO_KEYS and COMMON[k_] is None and not any(s.startswith(k_ + '=') for s in a.set)
            print(f'{k_:14s} {fmt(p[k_])}' + ('   (auto)' if auto and p[k_] is not None else ''))
        return
    probs = check_layout(p)
    if probs: sys.exit('cannot build this rack:\n  - ' + '\n  - '.join(probs))

    rows = p['_rows']
    types = [R['type'] for R in rows]
    if a.out: out = a.out
    else:
        out = (f'column_rack_{types[0]}mL_{len(rows)}x{rows[0]["n"]}' if len(set(types)) == 1
               else 'column_rack_' + '-'.join(f'{R["type"]}x{R["n"]}' for R in rows)) \
              + ('_feet' if p['feet'] else '') + '.stl'

    U, F, T, vol = build_mesh(p)
    Ttxt, tinfo = text_triangles(p) if p['text'] else (np.zeros((0, 3, 3)), None)
    write_stl(out, np.concatenate([T, Ttxt]) if len(Ttxt) else T)
    if a.profile_csv:
        R = rows[0]
        xs = [(R['centres'][0] + R['centres'][1]) / 2 if R['n'] > 1 else R['centres'][0] + R['pitch'] / 2,
              R['centres'][0]]
        with open(a.profile_csv, 'w') as f:
            f.write('where,y,z\n')
            for name, x in zip(('finger', 'slot_centre'), xs):
                for y, z in section(p, x): f.write(f'{name},{y:.4f},{z:.4f}\n')
    if a.preview: preview(a.preview, np.concatenate([T, Ttxt]) if len(Ttxt) else T, p)

    W, H, D = p['_W'], p['_H'], p['shelf_depth'] + p['plate_t']
    mass = vol / 1000 * p['density']
    print(f'{out}: {len(T)} triangles, volume {vol/1000:.1f} cm3 (~{mass:.0f} g solid at {p["density"]} g/cm3)')
    print(f'  size {W:.1f} (W) x {D:.1f} (D) x {H:.1f} (H) mm; '
          f'closed manifold: {check_manifold(F)}, zero-area triangles: {sliver_count(T)}, '
          f'cross-sections simple: {sections_simple(p)}')
    for i, R in enumerate(rows):
        print(f'  row {i + 1} (top first): {R["n"]} x {R["type"]} mL, slot {R["w"]:g} mm, pitch {R["pitch"]:g} mm, '
              f'slot centre {R["yc"]:.1f} mm behind the front edge, shelf top z={R["zt"]:.1f}')
    if p['hook_n']:
        hc = hole_centres(p)
        shape = f'{p["hook_d"] + p["hook_slot"]:g} x {p["hook_d"]:g} mm slots' if p['hook_slot'] else f'Ø{p["hook_d"]:g} mm'
        print(f'  hook holes: {len(hc)} x {shape} at x = {", ".join(f"{c:.1f}" for c in hc)}'
              + (f' (spacing {hc[-1] - hc[0]:.1f} mm' + (f', ±{p["hook_slot"] / 2:g} mm tolerance)' if p['hook_slot'] else ')')
                 if len(hc) > 1 else '')
              + f', {p["hook_from_top"]:g} mm below the top')
    if tinfo:
        print(f'  text "{p["text"]}": {tinfo["size"]:g} mm{" (auto)" if tinfo["auto"] else ""}, '
              f'{tinfo["x"][1] - tinfo["x"][0]:.0f} x {tinfo["z"][1] - tinfo["z"][0]:.0f} mm, raised {p["text_h"]:g} mm '
              f'(separate shells, merged by the slicer)')
        for pr in check_text(tinfo): print('  WARNING: ' + pr)
    # printer bed (printed lying on its back)
    margin = max(min(p['bed_x'] - W, p['bed_y'] - H), min(p['bed_y'] - W, p['bed_x'] - H))
    fits = margin >= 0
    print(f'  print footprint {W:.0f} x {H:.0f} mm on its back: '
          + (f'fits a {p["bed_x"]:g} x {p["bed_y"]:g} mm bed ({margin:.0f} mm to spare'
             + (', tight: use no brim/skirt)' if margin < 5 else ')') if fits
             else f'WARNING: does not fit a {p["bed_x"]:g} x {p["bed_y"]:g} mm bed'))
    if not fits:
        R = rows[0]; nmax = int((max(p['bed_x'], p['bed_y']) - 2 * p['end_extra']) // R['pitch'])
        print(f'           (at most {nmax} slots per row at {R["pitch"]:g} mm pitch: use --width {max(p["bed_x"], p["bed_y"]):g}, '
              f'or set bed_x/bed_y for a larger printer)')
    # load
    ncol = sum(R['n'] for R in rows)
    if p['col_mass']:
        print(f'  full load: {ncol} columns x {p["col_mass"]:g} g = {ncol * p["col_mass"]:.0f} g '
              f'+ rack ~{mass:.0f} g (100 % infill) -> {(ncol * p["col_mass"] + mass) / 1000:.2f} kg'
              + (' (standing on its feet; if hung instead: on ' if p['feet'] else ' on ')
              + (f'{p["hook_n"]} hooks, {(ncol * p["col_mass"] + mass) / 1000 / max(p["hook_n"], 1):.2f} kg each'
                 if p['hook_n'] else 'no hooks')
              + (')' if p['feet'] else ''))
    # feet: standing against a wall
    if p['feet']:
        fs = foot_spans(p)
        print(f'  feet: {len(fs)} x {p["foot_w"]:g} mm, reaching {p["foot_ext"]:g} mm in front of the shelves '
              f'(footprint depth {D + p["foot_ext"]:.0f} mm); lowest shelf underside {rows[-1]["zb"]:.0f} mm above the table')
        # centre of mass (y) of the empty rack (solid) and of the full rack, vs the front of the feet
        Tm = T; vt = np.einsum('ij,ij->i', Tm[:, 0], np.cross(Tm[:, 1], Tm[:, 2])) / 6
        cy = (vt * Tm[:, :, 1].sum(1) / 4).sum() / vt.sum()
        msg = f'  centre of mass {cy + p["foot_ext"]:.0f} mm behind the front of the feet (empty)'
        if p['col_mass']:
            mc = sum(R['n'] for R in rows) * p['col_mass']; cyc = np.mean([R['yc'] for R in rows])
            cyf = (cy * mass + cyc * mc) / (mass + mc)
            msg += f', {cyf + p["foot_ext"]:.0f} mm (full)'
        print(msg + '; the wall stops it tipping backwards')
        free = rows[-1]['zb'] - p['foot_t']             # above a foot (feet may sit under a column end)
        if p['col_below'] and p['col_below'] > free - 2:
            print(f'  WARNING: bottom-row columns ({p["col_below"]:g} mm below the shelf) reach the table: '
                  f'increase foot_clear')
        if not p['col_below']:
            print(f'  note: bottom-row columns must hang less than {free - 2:.0f} mm below the shelf '
                  f'(give col_below, or set foot_clear)')
    # column clearances (optional)
    for i, R in enumerate(rows):
        tag = f'row {i + 1}'
        if p['shoulder_d']:
            if p['shoulder_d'] > R['pitch'] - 1:
                print(f'  WARNING {tag}: shoulder_d {p["shoulder_d"]:g} mm leaves < 1 mm between neighbouring columns (pitch {R["pitch"]:g})')
            room = p['shelf_depth'] - p['top_chamfer'] - R['yc']
            if p['shoulder_d'] / 2 > room:
                print(f'  WARNING {tag}: shoulder ({p["shoulder_d"] / 2:.1f} mm radius) hits the plate '
                      f'({room:.1f} mm from slot centre): increase back_gap')
            if p['shoulder_d'] <= R['w']:
                print(f'  WARNING {tag}: shoulder_d is not larger than the slot: the column would fall through')
        if i + 1 < len(rows) and p['col_above'] is not None and p['col_below'] is not None:
            lo = rows[i + 1]
            gap = (R['zb'] - p['col_below']) - (lo['zt'] + p['col_above'])
            need = p['lip_h'] + 2
            if gap < need:
                print(f'  WARNING rows {i + 1}/{i + 2}: only {gap:.1f} mm between a hanging column and the one below '
                      f'(need about {need:g} mm to lift a column over the lip): increase row_pitch')
            else:
                print(f'  rows {i + 1}/{i + 2}: {gap:.1f} mm free between columns')

if __name__ == '__main__':
    main()
