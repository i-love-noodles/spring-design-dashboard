import math

import numpy as np
import streamlit as st

# ── Material data (Shigley's, Table 10-4) ──
# Each entry: (A (psi), m, G (psi))
WIRE_MATERIALS = {
    "Music Wire (ASTM A228)":  (201_000, 0.145, 11_500_000),
    "Hard Drawn (ASTM A227)":  (140_000, 0.190, 11_500_000),
    "Stainless (ASTM A313)":   (169_000, 0.146, 10_000_000),
}
WIRE_MAT_NAMES = list(WIRE_MATERIALS.keys())

WIRE_SIZES = [0.072, 0.075, 0.08, 0.085, 0.091, 0.095, 0.1, 0.105, 0.11]

MM_PER_IN = 25.4
LBF_IN_TO_J = 0.1129848
FPS_PER_MPS = 3.28084

BLASTER_PRESETS = {
    "Custom": None,
    "BABP": (134.0, 39.0),
    "Hijinx": (175.0, 50.0),
    "Lynx": (137.0, 38.0),
    "Minx": (112.0, 38.0),
    "Lonx": (226.0, 76.0),
}
BLASTER_PRESET_NAMES = list(BLASTER_PRESETS.keys())

# Spring presets: (wire_d_mm, active_coils, length_mm, od_mm, end_type)
SPRING_PRESETS = {
    "Custom": None,
    "SFX3":   (2.4,  10.625, 145.5, 30.5, "Closed and ground"),
    "SFX4":   (2.3,  11.0,   144.0, 30.8, "Closed and ground"),
    "SFX5":   (2.1,  11.0,   146.5, 30.6, "Closed and ground"),
    "SFX6":   (2.0,  11.5,   141.5, 31.0, "Closed and ground"),
    "SFX7":   (1.9,  12.625, 143.0, 30.3, "Closed and ground"),
    "SFX8":   (1.9,  15.125, 144.5, 30.4, "Closed and ground"),
    "5kg LS": (2.0,  11.5,   145.0, 30.2, "Closed and ground"),
    "8kg LS": (2.3,  12.25,  141.0, 30.5, "Closed and ground"),
}
SPRING_PRESET_NAMES = list(SPRING_PRESETS.keys())


def qp(key, default):
    """Read a typed default from the URL query string, falling back to *default*."""
    val = st.query_params.get(key)
    if val is None:
        return default
    try:
        return type(default)(val)
    except (ValueError, TypeError):
        return default


def snap(val, lo, hi, step):
    """Clamp *val* to [lo, hi] and round to the nearest *step*."""
    val = max(lo, min(hi, val))
    return round(round((val - lo) / step) * step + lo, 10)


def linked(label, key, lo, hi, default, step, fmt="%.3f", container=None,
           slider_lo=None, slider_hi=None, help=None):
    """Render a linked number-input + slider pair, returning the current value."""
    ct = container or st
    s_lo = slider_lo if slider_lo is not None else lo
    s_hi = slider_hi if slider_hi is not None else hi
    nk, sk = f"{key}_n", f"{key}_s"

    def _sn():
        st.session_state[sk] = max(s_lo, min(s_hi, st.session_state[nk]))

    def _ns():
        st.session_state[nk] = st.session_state[sk]

    if nk not in st.session_state:
        st.session_state[nk] = default
        st.session_state[sk] = snap(default, s_lo, s_hi, step)

    c1, c2 = ct.columns([1, 2])
    with c1:
        st.number_input(
            label, min_value=lo, max_value=hi, step=step,
            format=fmt, key=nk, on_change=_sn,
            help=help,
        )
    with c2:
        st.slider(
            label, min_value=s_lo, max_value=s_hi, step=step,
            format=fmt, key=sk, on_change=_ns, label_visibility="collapsed",
        )
    return st.session_state[nk]


# ── Spring physics ──

def compute_spring(d, OD, Na, Lf, end_type, wire_type):
    """Compute all spring properties from basic inputs.  Returns a dict."""
    A, m, G = WIRE_MATERIALS[wire_type]
    D = OD - d
    C = D / d
    Kw = (4 * C - 1) / (4 * C - 4) + 0.615 / C
    k = (G * d**4) / (8 * D**3 * Na)

    Nt = Na + 2
    Hs = Nt * d if end_type == "Closed and ground" else (Nt + 1) * d

    Sut = A / d**m
    tau_allow = 0.45 * Sut
    F_safe = (tau_allow * math.pi * d**3) / (8 * D * Kw)
    x_safe = F_safe / k
    L_safe = Lf - x_safe
    L_safe_rpt = max(L_safe, Hs)
    safe_to_solid = L_safe < Hs

    x_solid = Lf - Hs
    F_solid = k * x_solid
    tau_solid = Kw * (8 * F_solid * D) / (math.pi * d**3)
    util_solid = tau_solid / tau_allow

    ID = OD - 2 * d
    pitch = (Lf - 2 * d) / Na if end_type == "Closed and ground" else (Lf - 3 * d) / Na
    coils_per_inch = 1.0 / pitch if pitch > 0 else 0.0
    sp_between = pitch - d
    weight_lbf = (math.pi * D * Nt) * (math.pi * (d / 2) ** 2) * 0.284

    return dict(
        d=d, OD=OD, Na=Na, Lf=Lf, end_type=end_type, wire_type=wire_type,
        A=A, m=m, G=G, D=D, C=C, Kw=Kw, k=k, Nt=Nt, Hs=Hs,
        Sut=Sut, tau_allow=tau_allow, F_safe=F_safe, x_safe=x_safe,
        L_safe=L_safe, L_safe_rpt=L_safe_rpt, safe_to_solid=safe_to_solid,
        x_solid=x_solid, F_solid=F_solid, tau_solid=tau_solid, util_solid=util_solid,
        ID=ID, pitch=pitch, coils_per_inch=coils_per_inch,
        sp_between=sp_between, weight_lbf=weight_lbf,
    )


def find_candidates(*, target_mode, target_fps=None, target_rate=None,
                    efficiency=0.50, dart_kg=0.0012,
                    comp_from_mm, comp_to_mm, margin_mm=2.0,
                    od_mode="fixed", od_fixed=1.4, od_min=1.0, od_max=2.0,
                    Lf, wire_type, end_type):
    """Search wire diameters and OD values for feasible spring designs.

    Returns (candidates, reject_reasons) where candidates is a sorted list of
    dicts and reject_reasons tallies why designs were rejected.
    """
    A, m, G = WIRE_MATERIALS[wire_type]
    comp_from_in = comp_from_mm / MM_PER_IN
    comp_to_in = comp_to_mm / MM_PER_IN

    if od_mode == "fixed":
        od_values = [od_fixed]
    else:
        od_values = list(np.arange(od_min, od_max + 0.005, 0.01))

    candidates = []
    rejects = {"spring_index": 0, "solid_too_tall": 0, "overstressed": 0, "na_too_low": 0}

    for d in WIRE_SIZES:
        Sut = A / d**m
        tau_allow = 0.45 * Sut

        for OD in od_values:
            OD = round(OD, 3)
            D = OD - d
            if D <= 0:
                continue
            C = D / d
            if C < 3 or C > 16:
                rejects["spring_index"] += 1
                continue

            Kw = (4 * C - 1) / (4 * C - 4) + 0.615 / C

            # Determine required spring rate
            if target_mode == "fps":
                v_mps = target_fps / FPS_PER_MPS
                E_needed = 0.5 * dart_kg * v_mps**2 / efficiency
                x_from = Lf - comp_from_in
                x_to = Lf - comp_to_in
                denom = x_to**2 - x_from**2
                if denom <= 0:
                    continue
                k_target = 2.0 * (E_needed / LBF_IN_TO_J) / denom
            else:
                k_target = target_rate

            if k_target <= 0:
                continue

            # Solve for Na, round to 0.25
            Na_exact = (G * d**4) / (8 * D**3 * k_target)
            Na = round(Na_exact * 4) / 4
            if Na < 1.0:
                rejects["na_too_low"] += 1
                continue

            k_actual = (G * d**4) / (8 * D**3 * Na)

            Nt = Na + 2
            Hs = Nt * d if end_type == "Closed and ground" else (Nt + 1) * d
            Hs_mm = Hs * MM_PER_IN

            if Hs_mm > comp_to_mm:
                rejects["solid_too_tall"] += 1
                continue

            # Stress at compress-to
            x_at_ct = Lf - comp_to_in
            F_at_ct = k_actual * x_at_ct
            tau_at_ct = Kw * (8 * F_at_ct * D) / (math.pi * d**3)
            util_at_ct = tau_at_ct / tau_allow

            if util_at_ct > 1.0:
                rejects["overstressed"] += 1
                continue

            # Stress at solid
            x_at_solid = Lf - Hs
            F_at_solid = k_actual * x_at_solid
            tau_at_solid = Kw * (8 * F_at_solid * D) / (math.pi * d**3)
            util_at_solid = tau_at_solid / tau_allow
            safe_to_solid = util_at_solid <= 1.0

            # Compute actual FPS
            x_from = Lf - comp_from_in
            x_to = Lf - comp_to_in
            E_stored_in = 0.5 * k_actual * (x_to**2 - x_from**2)
            E_stored_j = E_stored_in * LBF_IN_TO_J
            if dart_kg > 0 and E_stored_j > 0:
                fps_actual = math.sqrt(2 * efficiency * E_stored_j / dart_kg) * FPS_PER_MPS
            else:
                fps_actual = 0.0

            margin = comp_to_mm - Hs_mm
            index_score = abs(C - 8.5)

            candidates.append(dict(
                d=d, OD=OD, Na=Na, k=k_actual, Hs_in=Hs, Hs_mm=Hs_mm,
                margin_mm=margin, safe_to_solid=safe_to_solid,
                util_at_ct=util_at_ct, util_at_solid=util_at_solid,
                fps=fps_actual, Lf=Lf, C=C, Nt=Nt, Kw=Kw,
                index_score=index_score,
            ))

    # Rank: safe-to-solid first, then largest margin, then best spring index
    candidates.sort(key=lambda c: (not c["safe_to_solid"], -c["margin_mm"], c["index_score"]))

    return candidates, rejects


def pareto_filter(candidates):
    """Return only Pareto-optimal candidates on (margin_mm max, util_at_ct min)."""
    result = []
    for c in candidates:
        dominated = any(
            o["margin_mm"] >= c["margin_mm"] and o["util_at_ct"] <= c["util_at_ct"]
            and (o["margin_mm"] > c["margin_mm"] or o["util_at_ct"] < c["util_at_ct"])
            for o in candidates if o is not c
        )
        if not dominated:
            result.append(c)
    return result
