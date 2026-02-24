import math

import streamlit as st

from spring_helpers import (
    WIRE_MATERIALS, WIRE_MAT_NAMES, WIRE_SIZES_IN, WIRE_SIZES_MM,
    MM_PER_IN, LBF_IN_TO_J, FPS_PER_MPS,
    BLASTER_PRESETS, BLASTER_PRESET_NAMES,
    SPRING_PRESETS, SPRING_PRESET_NAMES,
    END_TYPE_NAMES,
    compute_spring, format_wire_d, linked, qp, best_blaster_for_spring,
)

st.title("Compression Spring Design Analyzer")

# ── Spring parameter inputs ──

st.markdown(
    '<div style="display:flex;align-items:center;gap:1.5rem;flex-wrap:wrap;margin-bottom:0.5rem">'
    '<h3 style="margin:0">Spring Parameters</h3>'
    '</div>',
    unsafe_allow_html=True,
)
_ctrl_cols = st.columns([1, 2, 1, 1, 1])
_ctrl_cols[0].markdown("**Units:**")
_unit_systems = ["Imperial", "Metric"]
_unit_sys_def = qp("units", "Imperial")
_unit_sys_idx = _unit_systems.index(_unit_sys_def) if _unit_sys_def in _unit_systems else 0
unit_system = _ctrl_cols[1].radio("Units", _unit_systems, index=_unit_sys_idx,
                                  horizontal=True, label_visibility="collapsed", key="unit_system")
_metric = unit_system == "Metric"

_prev_unit = st.session_state.get("_prev_unit")
if _prev_unit is not None and _prev_unit != unit_system:
    for _uk in ("od", "lf"):
        _nk, _sk = f"{_uk}_n", f"{_uk}_s"
        if _nk in st.session_state:
            if _metric:
                st.session_state[_nk] = st.session_state[_nk] * MM_PER_IN
                st.session_state[_sk] = st.session_state[_sk] * MM_PER_IN
            else:
                st.session_state[_nk] = st.session_state[_nk] / MM_PER_IN
                st.session_state[_sk] = st.session_state[_sk] / MM_PER_IN
st.session_state["_prev_unit"] = unit_system

_ctrl_cols[2].markdown("**Wire Sizes:**",
                       help="Which wire diameter standards to include in the dropdown. "
                            "Check both to see imperial and metric sizes interleaved by diameter.")

_wire_auto_msg = None
_pending_preset = st.session_state.get("spring_preset",
                                       st.query_params.get("preset", "Custom"))
if _pending_preset in SPRING_PRESETS and SPRING_PRESETS[_pending_preset] is not None:
    _pp = SPRING_PRESETS[_pending_preset]
    _pp_is_metric = any(abs(_pp[0] - mm) < 0.01 for mm in WIRE_SIZES_MM)
    if _pp_is_metric and not st.session_state.get("wire_inc_mm", False):
        st.session_state["wire_inc_mm"] = True
        _wire_auto_msg = "Enabled metric wire sizes to match preset."
    elif not _pp_is_metric and not st.session_state.get("wire_inc_in", True):
        st.session_state["wire_inc_in"] = True
        _wire_auto_msg = "Enabled imperial wire sizes to match preset."

_inc_in = _ctrl_cols[3].checkbox("in", value=True, key="wire_inc_in")
_inc_mm = _ctrl_cols[4].checkbox("mm", value=False, key="wire_inc_mm")

_preset_def = st.query_params.get("preset", "Custom")
_preset_idx = SPRING_PRESET_NAMES.index(_preset_def) if _preset_def in SPRING_PRESET_NAMES else 0
spring_preset = st.selectbox("Spring Preset", SPRING_PRESET_NAMES, index=_preset_idx,
                             key="spring_preset",
                             help="Select a known spring to auto-fill parameters.")
if _wire_auto_msg:
    st.info(_wire_auto_msg)
if spring_preset != "Custom":
    _sp = SPRING_PRESETS[spring_preset]
    _sp_d_mm = _sp[0]
    _sp_is_metric = any(abs(_sp_d_mm - mm) < 0.01 for mm in WIRE_SIZES_MM)
    if _sp_is_metric:
        _sp_d_in = _sp_d_mm / MM_PER_IN
    else:
        _sp_d_in = min(WIRE_SIZES_IN, key=lambda w: abs(w - _sp_d_mm / MM_PER_IN))
    if _metric:
        st.session_state["od_n"] = round(_sp[3], 1)
        st.session_state["od_s"] = round(_sp[3], 1)
        st.session_state["lf_n"] = round(_sp[2], 1)
        st.session_state["lf_s"] = round(_sp[2], 1)
    else:
        _sp_od_in = round(_sp[3] / MM_PER_IN, 4)
        _sp_lf_in = round(_sp[2] / MM_PER_IN, 4)
        st.session_state["od_n"] = _sp_od_in
        st.session_state["od_s"] = _sp_od_in
        st.session_state["lf_n"] = _sp_lf_in
        st.session_state["lf_s"] = _sp_lf_in
    st.session_state["na_n"] = _sp[1]
    st.session_state["na_s"] = _sp[1]
    _sp_end = _sp[4]

    if st.session_state.get("_prev_spring_preset") != spring_preset:
        _best_bl = best_blaster_for_spring(_sp[2])
        st.session_state["analysis_preset"] = _best_bl
        if _best_bl != "Custom":
            _bl_from, _bl_to = BLASTER_PRESETS[_best_bl]
            st.session_state["comp_from_n"] = _bl_from
            st.session_state["comp_from_s"] = _bl_from
            st.session_state["comp_to_n"] = _bl_to
            st.session_state["comp_to_s"] = _bl_to
st.session_state["_prev_spring_preset"] = spring_preset

p_left, p_mid, p_right = st.columns(3)

with p_left:
    if _metric:
        _od_mm = linked("Outer Diameter (mm)", "od", 5.0, 127.0,
                         qp("od", 1.40) * MM_PER_IN, 0.1, "%.1f", container=p_left,
                         help="Outside diameter of the coil.",
                         preset_key="spring_preset")
        OD = _od_mm / MM_PER_IN
        _lf_mm = linked("Free Length (mm)", "lf", 12.0, 508.0,
                         qp("lf", 5.60) * MM_PER_IN, 1.0, "%.1f", container=p_left,
                         help="Length of the spring with no load applied.",
                         preset_key="spring_preset")
        Lf = _lf_mm / MM_PER_IN
    else:
        OD = linked("Outer Diameter (in)", "od", 0.20, 5.00, qp("od", 1.40), 0.005,
                    container=p_left,
                    help="Outside diameter of the coil (not the mean diameter).",
                    preset_key="spring_preset")
        Lf = linked("Free Length (in)", "lf", 0.50, 20.00, qp("lf", 5.60), 0.05, "%.2f",
                    container=p_left,
                    help="Length of the spring with no load applied.",
                    preset_key="spring_preset")

# ── Wire size list from checkboxes ──

_mm_as_in = [mm / MM_PER_IN for mm in WIRE_SIZES_MM]
_wire_list = []
if _inc_in:
    _wire_list.extend(WIRE_SIZES_IN)
if _inc_mm:
    _wire_list.extend(_mm_as_in)
_wire_list = sorted(set(_wire_list))

if not _wire_list:
    st.error("Select at least one wire size group (in or mm).")
    st.stop()

with p_mid:
    if spring_preset != "Custom":
        _d_def = _sp_d_in
    else:
        _d_def = min(_wire_list, key=lambda w: abs(w - float(qp("d", 0.095))))

    _d_idx = min(range(len(_wire_list)), key=lambda i: abs(_wire_list[i] - _d_def))
    st.markdown("**Wire Diameter**")
    d = st.selectbox(
        "Wire Diameter", _wire_list, index=_d_idx,
        format_func=format_wire_d, label_visibility="collapsed",
        help="Cross-section diameter of the spring wire.",
    )
    Na = linked("Number of Active Coils", "na", 1.0, 50.0, qp("na", 10.0), 0.25,
                "%.2f", container=p_mid,
                help="Coils that deflect under load. Excludes the closed end coils.",
                preset_key="spring_preset")

with p_right:
    if spring_preset != "Custom":
        _end_def = _sp_end
    else:
        _end_def = st.query_params.get("end", END_TYPE_NAMES[0])
    _end_idx = END_TYPE_NAMES.index(_end_def) if _end_def in END_TYPE_NAMES else 0
    st.markdown("**End Type**")
    end_type = st.selectbox("End Type", END_TYPE_NAMES, index=_end_idx,
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
           "na": str(Na), "end": end_type, "preset": spring_preset,
           "units": unit_system}
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
    _PSI_TO_MPA = 0.00689476
    _LBF_TO_N = 4.44822
    _D_mm = D * MM_PER_IN
    _OD_mm = OD * MM_PER_IN
    _d_mm = d * MM_PER_IN
    _k_n_mm = k * _LBF_TO_N / MM_PER_IN
    _Sut_mpa = Sut * _PSI_TO_MPA
    _tau_allow_mpa = tau_allow * _PSI_TO_MPA
    _F_safe_n = F_safe * _LBF_TO_N
    _x_safe_mm = x_safe * MM_PER_IN
    _F_solid_n = F_solid * _LBF_TO_N
    _tau_solid_mpa = tau_solid * _PSI_TO_MPA
    _G_mpa = G * _PSI_TO_MPA
    _G_gpa = _G_mpa / 1000

    dc_left, dc_right = st.columns(2)
    with dc_left:
        st.markdown(
            "| Property | Value |\n|---|---|\n"
            f"| Mean Diameter (D) | {_D_mm:.2f} mm |\n"
            f"| Spring Index (C) | {C:.2f} |\n"
            f"| Wahl Factor (Kw) | {Kw:.4f} |\n"
            f"| Total Coils (Nt) | {Nt:g} |\n"
            f"| Sut (tensile) | {_Sut_mpa:.0f} MPa |\n"
            f"| &tau;_allow (45%) | {_tau_allow_mpa:.0f} MPa |\n"
            f"| Safe Load (F_safe) | {_F_safe_n:.2f} N |\n"
            f"| Safe Deflection | {_x_safe_mm:.1f} mm |\n"
            f"| Shear Stress at Solid | {_tau_solid_mpa:.0f} MPa |\n"
            f"| Shear Modulus (G) | {_G_gpa:.1f} GPa |\n"
        )
    with dc_right:
        st.markdown("**Formulas**")
        st.latex(rf"D = OD - d = {_OD_mm:.2f} - {_d_mm:.2f} = {_D_mm:.2f}\;\text{{mm}}")
        st.latex(rf"C = D / d = {_D_mm:.2f} / {_d_mm:.2f} = {C:.2f}")
        st.latex(rf"K_w = \frac{{4C-1}}{{4C-4}} + \frac{{0.615}}{{C}} = {Kw:.4f}")
        st.latex(rf"k = \frac{{G \, d^4}}{{8 \, D^3 \, N_a}} = \frac{{{_G_mpa:.0f} \times {_d_mm:.2f}^4}}{{8 \times {_D_mm:.2f}^3 \times {Na:g}}} = {_k_n_mm:.3f}\;\text{{N/mm}}")
        st.latex(rf"S_{{ut}} = \frac{{A}}{{d^m}} = {_Sut_mpa:.0f}\;\text{{MPa}}")
        st.latex(rf"\tau_{{allow}} = 0.45 \, S_{{ut}} = {_tau_allow_mpa:.0f}\;\text{{MPa}}")
        st.latex(rf"F_{{safe}} = \frac{{\tau_{{allow}} \, \pi \, d^3}}{{8 \, D \, K_w}} = {_F_safe_n:.2f}\;\text{{N}}")
        st.latex(rf"x_{{safe}} = F_{{safe}} / k = {_x_safe_mm:.1f}\;\text{{mm}}")
        st.latex(rf"\tau_{{solid}} = K_w \frac{{8 \, F_{{solid}} \, D}}{{\pi \, d^3}} = {_tau_solid_mpa:.0f}\;\text{{MPa}}")

# ── Operating-point check (optional) ──

if st.checkbox("Show Operating-Point Check"):
    _LBF_TO_N = 4.44822
    if _metric:
        _Hs_mm = round(Hs * MM_PER_IN, 1)
        _Lf_mm = round(Lf * MM_PER_IN, 1)
        _def_mm = round(max(_Lf_mm * 0.7, _Hs_mm + 0.1), 1)
        L_check_mm = st.slider(
            "Loaded length (mm)", min_value=_Hs_mm, max_value=_Lf_mm,
            value=min(_def_mm, _Lf_mm), step=0.1, format="%.1f",
        )
        L_check = L_check_mm / MM_PER_IN
    else:
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
    if _metric:
        c1.metric("Deflection", f"{x_op * MM_PER_IN:.1f} mm")
        c2.metric("Force", f"{F_op * _LBF_TO_N:.2f} N")
        c3.metric("Shear Stress", f"{tau_op / 145.038:.0f} MPa")
    else:
        c1.metric("Deflection", f"{x_op:.3f} in")
        c2.metric("Force", f"{F_op:.2f} lbf")
        c3.metric("Shear Stress", f"{tau_op:,.0f} psi")
    c4.metric(
        "Utilization", f"{util:.1%}",
        delta="SAFE" if util <= 1 else "OVER",
        delta_color="normal" if util <= 1 else "inverse",
    )

# ── Blaster estimation ──

_url_comp = st.query_params.get("comp_from")
if _url_comp is not None:
    _opt_sig = f"{_url_comp}_{st.query_params.get('comp_to')}_{st.query_params.get('blaster')}"
    if st.session_state.get("_opt_sig") != _opt_sig:
        st.session_state["_opt_sig"] = _opt_sig
        _url_bl = st.query_params.get("blaster", "Custom")
        if _url_bl in BLASTER_PRESET_NAMES:
            st.session_state["analysis_preset"] = _url_bl
        st.session_state["comp_from_n"] = qp("comp_from", 134.0)
        st.session_state["comp_from_s"] = qp("comp_from", 134.0)
        st.session_state["comp_to_n"] = qp("comp_to", 39.0)
        st.session_state["comp_to_s"] = qp("comp_to", 39.0)

show_fps = st.checkbox("Estimate FPS", value=True)
show_eff = st.checkbox("Estimate Efficiency")

if show_fps or show_eff:
    bl_ct = st.container()

    _bl_def = st.query_params.get("blaster", "Lynx")
    _bl_idx = BLASTER_PRESET_NAMES.index(_bl_def) if _bl_def in BLASTER_PRESET_NAMES else 0
    preset = bl_ct.selectbox("Blaster Preset", BLASTER_PRESET_NAMES, index=_bl_idx,
                             key="analysis_preset",
                             help="Select a blaster to auto-fill compress from/to values.")
    bl_ct.caption("This checks spring performance only — it does not verify physical "
                  "fit with plunger tube ID, barrel length, or other blaster geometry.")
    if preset != "Custom":
        _preset_from, _preset_to = BLASTER_PRESETS[preset]
        st.session_state["comp_from_n"] = _preset_from
        st.session_state["comp_from_s"] = _preset_from
        st.session_state["comp_to_n"] = _preset_to
        st.session_state["comp_to_s"] = _preset_to

    bl_left, bl_right = bl_ct.columns(2)

    dart_g = linked(
        "Dart Weight (g)", "dart_g", 0.50, 10.00, 1.00, 0.01, "%.2f",
        container=bl_left, slider_lo=0.80, slider_hi=1.30,
        help="Mass of the dart in grams.",
    )
    dart_kg = dart_g / 1000.0

    comp_from = linked(
        "Compression From (mm)", "comp_from", 5.0, 500.0,
        qp("comp_from", 134.0), 0.5, "%.1f",
        container=bl_left,
        help="Spring length when the blaster is unprimed (plunger forward, at rest).",
        preset_key="analysis_preset",
    )
    comp_to = linked(
        "Compression To (mm)", "comp_to", 5.0, 500.0,
        qp("comp_to", 39.0), 0.5, "%.1f",
        container=bl_right,
        help="Spring length when the blaster is primed (plunger pulled back, spring fully compressed).",
        preset_key="analysis_preset",
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
        if Lf * MM_PER_IN > comp_from:
            _preload_mm = Lf * MM_PER_IN - comp_from
            _preload_force = k * x_from
            st.info(
                f"Spring free length ({Lf * MM_PER_IN:.1f} mm) exceeds compress-from "
                f"({comp_from:.1f} mm) — the spring is preloaded by "
                f"{_preload_mm:.1f} mm ({_preload_force:.2f} lbf) at rest."
            )

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
            st.caption("FPS estimates are very rough approximations and should be treated "
                       "as a vibe check, not precise predictions.")

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
