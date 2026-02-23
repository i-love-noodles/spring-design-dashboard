import math

import streamlit as st

from spring_helpers import (
    WIRE_MATERIALS, WIRE_MAT_NAMES, WIRE_SIZES, MM_PER_IN, LBF_IN_TO_J, FPS_PER_MPS,
    BLASTER_PRESETS, BLASTER_PRESET_NAMES,
    SPRING_PRESETS, SPRING_PRESET_NAMES,
    compute_spring, linked, qp,
)

st.title("Compression Spring Design Analyzer")

# ── Spring parameter inputs ──

st.subheader("Spring Parameters")

_preset_def = st.query_params.get("preset", "Custom")
_preset_idx = SPRING_PRESET_NAMES.index(_preset_def) if _preset_def in SPRING_PRESET_NAMES else 0
spring_preset = st.selectbox("Spring Preset", SPRING_PRESET_NAMES, index=_preset_idx,
                             key="spring_preset",
                             help="Select a known spring to auto-fill parameters.")
if spring_preset != "Custom":
    _sp = SPRING_PRESETS[spring_preset]
    _sp_d_in = min(WIRE_SIZES, key=lambda w: abs(w - _sp[0] / MM_PER_IN))
    _sp_od_in = round(_sp[3] / MM_PER_IN, 4)
    _sp_lf_in = round(_sp[2] / MM_PER_IN, 4)
    st.session_state["od_n"] = _sp_od_in
    st.session_state["od_s"] = _sp_od_in
    st.session_state["lf_n"] = _sp_lf_in
    st.session_state["lf_s"] = _sp_lf_in
    st.session_state["na_n"] = _sp[1]
    st.session_state["na_s"] = _sp[1]

p_left, p_mid, p_right = st.columns(3)

with p_left:
    OD = linked("Outer Diameter (in)", "od", 0.20, 5.00, qp("od", 1.40), 0.005,
                container=p_left,
                help="Outside diameter of the coil (not the mean diameter).")
    Lf = linked("Free Length (in)", "lf", 0.50, 20.00, qp("lf", 5.60), 0.05, "%.2f",
                container=p_left,
                help="Length of the spring with no load applied.")

with p_mid:
    if spring_preset != "Custom":
        _d_def = _sp_d_in
    else:
        _d_def = qp("d", 0.095)
    _d_idx = WIRE_SIZES.index(_d_def) if _d_def in WIRE_SIZES else WIRE_SIZES.index(0.095)
    st.markdown("**Wire Diameter (in)**")
    d = st.selectbox(
        "Wire Diameter", WIRE_SIZES, index=_d_idx,
        format_func=lambda x: f"{x:.3f}", label_visibility="collapsed",
        help="Cross-section diameter of the spring wire.",
    )
    Na = linked("Number of Active Coils", "na", 1.0, 40.0, qp("na", 10.0), 0.25,
                "%.2f", container=p_mid,
                help="Coils that deflect under load. Excludes the closed end coils.")

with p_right:
    END_TYPES = ["Closed and ground", "Closed not ground"]
    _end_def = st.query_params.get("end", END_TYPES[0])
    _end_idx = END_TYPES.index(_end_def) if _end_def in END_TYPES else 0
    st.markdown("**End Type**")
    end_type = st.selectbox("End Type", END_TYPES, index=_end_idx,
                            label_visibility="collapsed",
                            help="'Closed and ground' ends are flat and squared off. "
                                 "'Closed not ground' ends are closed but not machined flat.")

    _mat_def = st.query_params.get("mat", WIRE_MAT_NAMES[0])
    _mat_idx = WIRE_MAT_NAMES.index(_mat_def) if _mat_def in WIRE_MAT_NAMES else 0
    st.markdown("**Wire Type**")
    wire_type = st.selectbox("Wire Type", WIRE_MAT_NAMES, index=_mat_idx,
                             label_visibility="collapsed",
                             help="Wire material grade. Affects tensile strength (Sut) and shear modulus (G).")

# ── Sync current values back to URL ──

_new_qp = {"od": f"{OD:.3f}", "lf": f"{Lf:.2f}", "mat": wire_type, "d": str(d),
           "na": str(Na), "end": end_type, "preset": spring_preset}
if any(st.query_params.get(k) != v for k, v in _new_qp.items()):
    st.query_params.update(**_new_qp)

# ── Core calculations ──

_s = compute_spring(d, OD, Na, Lf, end_type, wire_type)
_A, _M, G = _s["A"], _s["m"], _s["G"]
D, C, Kw, k = _s["D"], _s["C"], _s["Kw"], _s["k"]
Nt, Hs = _s["Nt"], _s["Hs"]
Sut, tau_allow = _s["Sut"], _s["tau_allow"]
F_safe, x_safe, L_safe, L_safe_rpt = _s["F_safe"], _s["x_safe"], _s["L_safe"], _s["L_safe_rpt"]
safe_to_solid = _s["safe_to_solid"]
x_solid, F_solid, tau_solid, util_solid = _s["x_solid"], _s["F_solid"], _s["tau_solid"], _s["util_solid"]
ID, pitch, coils_per_inch = _s["ID"], _s["pitch"], _s["coils_per_inch"]
sp_between, weight_lbf = _s["sp_between"], _s["weight_lbf"]

if D <= 0:
    st.error("Wire diameter >= outer diameter — invalid geometry.")
    st.stop()
if C < 3:
    st.warning(f"Spring index C = {C:.1f} is very low (< 3). Manufacturing may be difficult.")
elif C > 16:
    st.warning(f"Spring index C = {C:.1f} is very high (> 16). Spring may buckle or tangle.")

st.divider()

# ── Report ──

st.subheader("Calculated Specs")
x_safe_defl = Lf - L_safe_rpt
F_safe_at_rpt = k * x_safe_defl
spec_l, spec_m, spec_r = st.columns(3)

with spec_l:
    st.markdown(
        "| Spec | Value |\n|---|---|\n"
        f"| Spring Rate | {k:.3f} lbf/in |\n"
        f"| Inner Diameter | {ID:.3f} in |\n"
        f"| Total Coils | {Nt:g} |\n"
        f"| Spring Index | {C:.2f} |\n"
    )

with spec_m:
    st.markdown(
        "| Spec | Value |\n|---|---|\n"
        f"| Pitch | {pitch:.3f} in |\n"
        f"| Coils/Inch | {coils_per_inch:.3f} |\n"
        f"| Sp. Bet. Coils | {sp_between:.3f} in |\n"
        f"| Weight/Ea. | {weight_lbf:.4f} lbs |\n"
    )

with spec_r:
    _LBF_TO_KGF = 0.453592
    st.markdown(
        "| | Length | Load | Defl. |\n|---|---|---|---|\n"
        f"| **Safe** | {L_safe_rpt:.3f} in | {F_safe_at_rpt:.2f} lb | {x_safe_defl:.3f} in |\n"
        f"| | {L_safe_rpt * MM_PER_IN:.1f} mm | {F_safe_at_rpt * _LBF_TO_KGF:.2f} kg | {x_safe_defl * MM_PER_IN:.1f} mm |\n"
        f"| **Solid** | {Hs:.3f} in | {F_solid:.2f} lb | {x_solid:.3f} in |\n"
        f"| | {Hs * MM_PER_IN:.1f} mm | {F_solid * _LBF_TO_KGF:.2f} kg | {x_solid * MM_PER_IN:.1f} mm |\n"
    )

if safe_to_solid:
    st.info(f"Spring is safe to solid by the 45% rule — max {tau_solid / Sut:.0%} at solid.")
elif util_solid > 1.0:
    st.warning(f"Spring is overstressed at solid ({util_solid:.0%} of safe utilization) and may take a set.")
st.caption(f"\\*WB Jones 45% safe rule — {wire_type}, G = {G / 1e6:.1f} Mpsi")

# ── Detailed intermediate values ──

with st.expander("Detailed Calculations"):
    dc_left, dc_right = st.columns(2)
    with dc_left:
        st.markdown(
            "| Property | Value |\n|---|---|\n"
            f"| Mean Diameter (D) | {D:.4f} in |\n"
            f"| Spring Index (C) | {C:.2f} |\n"
            f"| Wahl Factor (Kw) | {Kw:.4f} |\n"
            f"| Total Coils (Nt) | {Nt:g} |\n"
            f"| Sut (tensile) | {Sut:,.0f} psi &ensp;({Sut / 1000:.1f} ksi) |\n"
            f"| &tau;_allow (45%) | {tau_allow:,.0f} psi &ensp;({tau_allow / 1000:.1f} ksi) |\n"
            f"| Safe Load (F_safe) | {F_safe:.3f} lbf |\n"
            f"| Safe Deflection | {x_safe:.3f} in |\n"
            f"| Shear Stress at Solid | {tau_solid:,.0f} psi &ensp;({tau_solid / 1000:.1f} ksi) |\n"
            f"| Shear Modulus (G) | {G / 1e6:.1f} Mpsi |\n"
        )
    with dc_right:
        st.markdown("**Formulas**")
        st.latex(rf"D = OD - d = {OD:.3f} - {d:.3f} = {D:.4f}")
        st.latex(rf"C = D / d = {D:.4f} / {d:.3f} = {C:.2f}")
        st.latex(rf"K_w = \frac{{4C-1}}{{4C-4}} + \frac{{0.615}}{{C}} = {Kw:.4f}")
        st.latex(rf"k = \frac{{G \, d^4}}{{8 \, D^3 \, N_a}} = \frac{{{G/1e6:.1f}\text{{M}} \times {d}^4}}{{8 \times {D:.4f}^3 \times {Na:g}}} = {k:.3f}\;\text{{lbf/in}}")
        st.latex(rf"S_{{ut}} = \frac{{A}}{{d^m}} = \frac{{{_A:,}}}{{{d}^{{{_M}}}}} = {Sut:,.0f}\;\text{{psi}}")
        st.latex(rf"\tau_{{allow}} = 0.45 \, S_{{ut}} = {tau_allow:,.0f}\;\text{{psi}}")
        st.latex(rf"F_{{safe}} = \frac{{\tau_{{allow}} \, \pi \, d^3}}{{8 \, D \, K_w}} = {F_safe:.3f}\;\text{{lbf}}")
        st.latex(rf"x_{{safe}} = F_{{safe}} / k = {x_safe:.3f}\;\text{{in}}")
        st.latex(rf"\tau_{{solid}} = K_w \frac{{8 \, F_{{solid}} \, D}}{{\pi \, d^3}} = {tau_solid:,.0f}\;\text{{psi}}")

# ── Operating-point check (optional) ──

if st.checkbox("Show Operating-Point Check"):
    default_len = round(max(Lf * 0.7, Hs + 0.01), 2)
    L_check = st.slider(
        "Loaded length (in)", min_value=round(Hs, 3), max_value=round(Lf, 3),
        value=min(default_len, round(Lf, 3)), step=0.01, format="%.3f",
    )

    x_op = Lf - L_check
    F_op = k * x_op
    tau_op = Kw * (8 * F_op * D) / (math.pi * d**3)
    util = tau_op / tau_allow

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Deflection", f"{x_op:.3f} in")
    c2.metric("Force", f"{F_op:.2f} lbf")
    c3.metric("Shear Stress", f"{tau_op:,.0f} psi")
    c4.metric(
        "Utilization", f"{util:.1%}",
        delta="SAFE" if util <= 1 else "OVER",
        delta_color="normal" if util <= 1 else "inverse",
    )

# ── Blaster estimation ──


show_fps = st.checkbox("Estimate FPS")
show_eff = st.checkbox("Estimate Efficiency")

if show_fps or show_eff:
    bl_ct = st.container()

    preset = bl_ct.selectbox("Blaster Preset", BLASTER_PRESET_NAMES, key="analysis_preset",
                             help="Select a blaster to auto-fill compress from/to values.")
    if preset != "Custom":
        _preset_from, _preset_to = BLASTER_PRESETS[preset]
        st.session_state["comp_from_n"] = _preset_from
        st.session_state["comp_from_s"] = _preset_from
        st.session_state["comp_to_n"] = _preset_to
        st.session_state["comp_to_s"] = _preset_to

    bl_left, bl_right = bl_ct.columns(2)

    dart_g = linked(
        "Dart Weight (g)", "dart_g", 0.50, 10.00, 1.20, 0.01, "%.2f",
        container=bl_left, slider_lo=0.80, slider_hi=1.30,
        help="Mass of the dart in grams.",
    )
    dart_kg = dart_g / 1000.0

    comp_from = linked(
        "Compression From (mm)", "comp_from", 5.0, 500.0, 134.0, 0.5, "%.1f",
        container=bl_left,
        help="Spring length when the blaster is unprimed (plunger forward, at rest).",
    )
    comp_to = linked(
        "Compression To (mm)", "comp_to", 5.0, 500.0, 39.0, 0.5, "%.1f",
        container=bl_right,
        help="Spring length when the blaster is primed (plunger pulled back, spring fully compressed).",
    )

    x_from = (Lf * MM_PER_IN - comp_from) / MM_PER_IN
    x_to = (Lf * MM_PER_IN - comp_to) / MM_PER_IN

    Hs_mm = Hs * MM_PER_IN
    if x_from < 0 or x_to < 0:
        st.error("Compression lengths exceed free length — invalid configuration.")
    elif comp_from <= comp_to:
        st.error("'Compression From' must be longer than 'Compression To'.")
    elif comp_to < Hs_mm:
        st.error(f"'Compression To' ({comp_to:.1f} mm) is below solid height ({Hs_mm:.1f} mm) — spring cannot compress that far.")
    else:
        energy_in = 0.5 * k * (x_to**2 - x_from**2)
        energy_j = energy_in * LBF_IN_TO_J

        if show_fps:
            st.markdown("### FPS Estimation")
            fps_left, fps_right = st.columns(2)
            eff_min = linked(
                "Efficiency Min (%)", "eff_min", 1.0, 100.0, 45.0, 1.0, "%.0f",
                container=fps_left,
                help="Low estimate of how much spring energy transfers to dart kinetic energy.",
            ) / 100.0
            eff_max = linked(
                "Efficiency Max (%)", "eff_max", 1.0, 100.0, 55.0, 1.0, "%.0f",
                container=fps_right,
                help="High estimate of how much spring energy transfers to dart kinetic energy.",
            ) / 100.0

            fps_lo = math.sqrt(2 * eff_min * energy_j / dart_kg) * FPS_PER_MPS
            fps_hi = math.sqrt(2 * eff_max * energy_j / dart_kg) * FPS_PER_MPS

            c1, c2, c3 = st.columns(3)
            c1.metric("Stored Energy", f"{energy_j:.3f} J")
            c2.metric("FPS (low est.)", f"{fps_lo:.0f} fps")
            c3.metric("FPS (high est.)", f"{fps_hi:.0f} fps")

        if show_eff:
            st.markdown("### Efficiency Estimation")
            measured_fps = st.number_input(
                "Measured FPS", min_value=0.0, max_value=2000.0,
                value=150.0, step=1.0, format="%.0f",
                help="Chronograph reading to back-calculate blaster efficiency.",
            )
            v_mps = measured_fps / FPS_PER_MPS
            ke_j = 0.5 * dart_kg * v_mps**2
            eff_est = ke_j / energy_j if energy_j > 0 else 0.0

            c1, c2, c3 = st.columns(3)
            c1.metric("Stored Energy", f"{energy_j:.3f} J")
            c2.metric("Dart KE", f"{ke_j:.3f} J")
            c3.metric("Efficiency", f"{eff_est:.1%}")
