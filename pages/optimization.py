import pandas as pd
import streamlit as st

from spring_helpers import (
    WIRE_MATERIALS, WIRE_MAT_NAMES, WIRE_SIZES_IN, WIRE_SIZES_MM,
    MM_PER_IN, LBF_IN_TO_J, FPS_PER_MPS,
    BLASTER_PRESETS, BLASTER_PRESET_NAMES,
    END_TYPE_NAMES,
    find_candidates, format_wire_d, pareto_filter, linked, qp,
)

st.title("Spring Design Optimization")

_unit_systems = ["Imperial", "Metric"]
_unit_sys_def = qp("units", "Imperial")
_unit_sys_idx = _unit_systems.index(_unit_sys_def) if _unit_sys_def in _unit_systems else 0
unit_system = st.radio("Units", _unit_systems, index=_unit_sys_idx,
                       horizontal=True, key="opt_unit_system")
_metric = unit_system == "Metric"

_prev_unit = st.session_state.get("_opt_prev_unit")
if _prev_unit is not None and _prev_unit != unit_system:
    for _uk in ("opt_od_fix", "opt_od_min", "opt_od_max"):
        _nk, _sk = f"{_uk}_n", f"{_uk}_s"
        if _nk in st.session_state:
            if _metric:
                st.session_state[_nk] = st.session_state[_nk] * MM_PER_IN
                st.session_state[_sk] = st.session_state[_sk] * MM_PER_IN
            else:
                st.session_state[_nk] = st.session_state[_nk] / MM_PER_IN
                st.session_state[_sk] = st.session_state[_sk] / MM_PER_IN
st.session_state["_opt_prev_unit"] = unit_system

# ══════════════════════════════════════════════
#  Inputs
# ══════════════════════════════════════════════

st.subheader("Target")
_tgt_modes = ["Target FPS", "Target Spring Rate"]
_tgt_def = st.query_params.get("tgt", _tgt_modes[0])
_tgt_idx = _tgt_modes.index(_tgt_def) if _tgt_def in _tgt_modes else 0
target_mode = st.radio("Optimize for", _tgt_modes, index=_tgt_idx,
                       horizontal=True, label_visibility="collapsed")

tgt_left, tgt_right = st.columns(2)

if target_mode == "Target FPS":
    with tgt_left:
        target_fps = linked("Target FPS", "opt_fps", 50.0, 500.0, qp("fps", 150.0), 1.0, "%.0f",
                            container=tgt_left,
                            help="Desired dart muzzle velocity in feet per second.")
        dart_g = linked("Dart Weight (g)", "opt_dart", 0.50, 10.00, qp("dart", 1.20), 0.01, "%.2f",
                        container=tgt_left, slider_lo=0.80, slider_hi=1.30,
                        help="Mass of the dart in grams.")
    with tgt_right:
        efficiency = linked("Efficiency (%)", "opt_eff", 1.0, 100.0, qp("eff", 50.0), 1.0, "%.0f",
                            container=tgt_right,
                            help="Estimated fraction of spring energy transferred to dart kinetic energy.") / 100.0
    dart_kg = dart_g / 1000.0
    target_rate = None
    _implied_rate_placeholder = st.empty()
else:
    with tgt_left:
        target_rate = linked("Spring Rate (lbf/in)", "opt_rate", 0.1, 100.0, qp("rate", 5.0), 0.1,
                             "%.2f", container=tgt_left,
                             help="Desired force per inch of deflection.")
    target_fps = None
    efficiency = 0.50
    dart_g = 1.20
    dart_kg = dart_g / 1000.0
    _implied_rate_placeholder = None

st.divider()

# ── Blaster geometry ──

st.subheader("Blaster Constraints")

_bp_def = st.query_params.get("blaster", "Custom")
_bp_idx = BLASTER_PRESET_NAMES.index(_bp_def) if _bp_def in BLASTER_PRESET_NAMES else 0
preset = st.selectbox("Blaster Preset", BLASTER_PRESET_NAMES, index=_bp_idx,
                      key="opt_preset",
                      help="Select a blaster to auto-fill compress from/to values.")

if preset != "Custom":
    _preset_from, _preset_to = BLASTER_PRESETS[preset]
    st.session_state["opt_cfrom_n"] = _preset_from
    st.session_state["opt_cfrom_s"] = _preset_from
    st.session_state["opt_cto_n"] = _preset_to
    st.session_state["opt_cto_s"] = _preset_to

bl_left, bl_mid, bl_right = st.columns(3)

comp_from_mm = linked("Compress From (mm)", "opt_cfrom", 5.0, 500.0, qp("cfrom", 134.0), 0.5, "%.1f",
                      container=bl_left,
                      help="Spring length when the blaster is unprimed (plunger forward, at rest).",
                      preset_key="opt_preset")
comp_to_mm = linked("Compress To (mm)", "opt_cto", 5.0, 500.0, qp("cto", 39.0), 0.5, "%.1f",
                    container=bl_mid,
                    help="Spring length when the blaster is primed (plunger pulled back, spring fully compressed).",
                    preset_key="opt_preset")
margin_mm = linked("Solid Height Margin (mm)", "opt_margin", 0.0, 20.0, qp("margin", 2.0), 0.5, "%.1f",
                   container=bl_right,
                   help="Extra clearance between the spring's solid height and compress-to. "
                        "Prevents coil bind and adds safety margin.")

st.divider()

# ── OD constraint ──

st.subheader("Spring Constraints")
c_left, c_mid, c_right = st.columns(3)

with c_left:
    _od_modes = ["Fixed OD", "OD Range"]
    _od_mode_def = st.query_params.get("od_mode", _od_modes[0])
    _od_mode_idx = _od_modes.index(_od_mode_def) if _od_mode_def in _od_modes else 0
    od_mode = st.radio("OD Constraint", _od_modes, index=_od_mode_idx, horizontal=True,
                       label_visibility="collapsed")
    if od_mode == "Fixed OD":
        if _metric:
            _od_fix_mm = linked("OD (mm)", "opt_od_fix", 5.0, 127.0,
                                qp("od_fix", 1.40) * MM_PER_IN, 0.1, "%.1f", container=c_left,
                                help="Outer diameter of the spring. Must fit inside the blaster plunger tube ID.")
            od_fixed = _od_fix_mm / MM_PER_IN
        else:
            od_fixed = linked("OD (in)", "opt_od_fix", 0.20, 5.00, qp("od_fix", 1.40), 0.005, container=c_left,
                              help="Outer diameter of the spring. Must fit inside the blaster plunger tube ID.")
        od_min = od_max = od_fixed
    else:
        if _metric:
            _od_min_mm = linked("OD Min (mm)", "opt_od_min", 5.0, 127.0,
                                qp("od_min", 1.20) * MM_PER_IN, 0.1, "%.1f", container=c_left,
                                help="Minimum outer diameter to search.")
            _od_max_mm = linked("OD Max (mm)", "opt_od_max", 5.0, 127.0,
                                qp("od_max", 1.50) * MM_PER_IN, 0.1, "%.1f", container=c_left,
                                help="Maximum outer diameter to search. Limited by blaster plunger tube ID.")
            od_min = _od_min_mm / MM_PER_IN
            od_max = _od_max_mm / MM_PER_IN
        else:
            od_min = linked("OD Min (in)", "opt_od_min", 0.20, 5.00, qp("od_min", 1.20), 0.01, container=c_left,
                            help="Minimum outer diameter to search.")
            od_max = linked("OD Max (in)", "opt_od_max", 0.20, 5.00, qp("od_max", 1.50), 0.01, container=c_left,
                            help="Maximum outer diameter to search. Limited by blaster plunger tube ID.")

# ── Free length ──

with c_mid:
    _lf_modes = ["From compress-from", "Exact"]
    _lf_mode_def = st.query_params.get("lf_mode", _lf_modes[0])
    _lf_mode_idx = _lf_modes.index(_lf_mode_def) if _lf_mode_def in _lf_modes else 0
    lf_mode = st.radio("Free Length", _lf_modes, index=_lf_mode_idx, horizontal=True,
                       label_visibility="collapsed")
    if lf_mode == "Exact":
        if _metric:
            _lf_mm = linked("Free Length (mm)", "opt_lf", 12.0, 508.0,
                            qp("opt_lf_val", 5.60) * MM_PER_IN, 1.0, "%.1f",
                            container=c_mid, slider_lo=comp_from_mm, slider_hi=280.0,
                            help="Length of the spring with no load applied.")
            Lf = _lf_mm / MM_PER_IN
        else:
            Lf = linked("Free Length (in)", "opt_lf", 0.50, 20.00, qp("opt_lf_val", 5.60), 0.05, "%.2f",
                        container=c_mid, slider_lo=comp_from_mm / MM_PER_IN, slider_hi=11.00,
                        help="Length of the spring with no load applied.")
    else:
        auto_settle = c_mid.checkbox("Auto estimate", value=True, key="opt_auto_settle",
                                     help="Automatically estimate settling allowance as ~5% of stroke. "
                                          "Only applies when the spring is safe — an overstressed spring "
                                          "will settle more and unpredictably.")
        if auto_settle:
            settling_mm = max(1.0, 0.05 * (comp_from_mm - comp_to_mm))
            c_mid.metric("Settling Allowance (mm)", f"{settling_mm:.1f}")
        else:
            settling_mm = linked("Settling Allowance (mm)", "opt_settle", 0.0, 25.0, 5.0, 0.5,
                                 "%.1f", container=c_mid,
                                 help="Extra length added to compress-from to account for the spring "
                                      "taking a set over time.")
        Lf = (comp_from_mm + settling_mm) / MM_PER_IN
        if _metric:
            c_mid.markdown(f"Free Length = **{Lf * MM_PER_IN:.1f} mm** ({Lf:.3f} in)")
        else:
            c_mid.markdown(f"Free Length = **{Lf:.3f} in** ({comp_from_mm + settling_mm:.1f} mm)")

# ── Wire type / end type ──

with c_right:
    st.markdown("**Wire Type**")
    wire_type = st.selectbox("Wire Type", WIRE_MAT_NAMES, label_visibility="collapsed",
                             key="opt_wire_type",
                             help="Wire material grade. Affects tensile strength (Sut) and shear modulus (G).")
    st.markdown("**End Type**")
    end_type = st.selectbox("End Type", END_TYPE_NAMES, label_visibility="collapsed",
                            key="opt_end_type",
                            help="'Closed and ground' ends are flat and squared off. "
                                 "'Closed not ground' ends are closed but not machined flat.")
    st.markdown("**Wire Sizes**")
    _inc_imperial = st.checkbox("Include imperial (in)", value=True, key="opt_inc_in")
    _inc_metric = st.checkbox("Include metric (mm)", value=False, key="opt_inc_mm")

_wire_sizes = []
if _inc_imperial:
    _wire_sizes.extend(WIRE_SIZES_IN)
if _inc_metric:
    _wire_sizes.extend([mm / MM_PER_IN for mm in WIRE_SIZES_MM])
_wire_sizes.sort()

if not _wire_sizes:
    st.error("Select at least one wire size group (imperial or metric).")
    st.stop()

# ── Sync current values back to URL ──

_new_qp = {"tgt": target_mode, "blaster": preset,
           "cfrom": f"{comp_from_mm:.1f}", "cto": f"{comp_to_mm:.1f}",
           "margin": f"{margin_mm:.1f}", "mat": wire_type, "end": end_type,
           "od_mode": od_mode, "lf_mode": lf_mode, "units": unit_system}
if target_mode == "Target FPS":
    _new_qp.update(fps=f"{target_fps:.0f}", dart=f"{dart_g:.2f}", eff=f"{efficiency * 100:.0f}")
else:
    _new_qp["rate"] = f"{target_rate:.2f}"
if od_mode == "Fixed OD":
    _new_qp["od_fix"] = f"{od_fixed:.3f}"
else:
    _new_qp.update(od_min=f"{od_min:.3f}", od_max=f"{od_max:.3f}")
if lf_mode == "Exact":
    _new_qp["opt_lf_val"] = f"{Lf:.2f}"

if any(st.query_params.get(k) != v for k, v in _new_qp.items()):
    st.query_params.update(**_new_qp)

# ══════════════════════════════════════════════
#  Validation
# ══════════════════════════════════════════════

if comp_from_mm <= comp_to_mm:
    st.error("'Compress From' must be longer than 'Compress To'.")
    st.stop()

if lf_mode == "Exact" and Lf < comp_from_mm / MM_PER_IN:
    st.error(f"Free length ({Lf:.2f} in / {Lf * MM_PER_IN:.1f} mm) must be longer than "
             f"compress-from ({comp_from_mm:.1f} mm). The spring needs to be longer than "
             f"its installed length.")
    st.stop()

if _implied_rate_placeholder is not None:
    v_mps = target_fps / FPS_PER_MPS
    E_needed_j = 0.5 * dart_kg * v_mps**2 / efficiency
    E_needed_in = E_needed_j / LBF_IN_TO_J
    x_from = Lf - comp_from_mm / MM_PER_IN
    x_to = Lf - comp_to_mm / MM_PER_IN
    denom = x_to**2 - x_from**2
    if denom > 0:
        k_implied = 2.0 * E_needed_in / denom
        _implied_rate_placeholder.info(
            f"Target FPS {target_fps:.0f} at {efficiency:.0%} efficiency "
            f"requires ~**{k_implied:.2f} lbf/in** spring rate"
        )

# ══════════════════════════════════════════════
#  Search
# ══════════════════════════════════════════════

st.divider()
st.subheader("Results")

with st.spinner("Searching for feasible designs..."):
    candidates, rejects = find_candidates(
        target_mode="fps" if target_mode == "Target FPS" else "rate",
        target_fps=target_fps,
        target_rate=target_rate,
        efficiency=efficiency,
        dart_kg=dart_kg,
        comp_from_mm=comp_from_mm,
        comp_to_mm=comp_to_mm,
        margin_mm=margin_mm,
        od_mode="fixed" if od_mode == "Fixed OD" else "range",
        od_fixed=od_fixed if od_mode == "Fixed OD" else None,
        od_min=od_min if od_mode != "Fixed OD" else None,
        od_max=od_max if od_mode != "Fixed OD" else None,
        Lf=Lf,
        wire_type=wire_type,
        end_type=end_type,
        wire_sizes=_wire_sizes,
    )

if not candidates:
    st.error("No feasible spring designs found.")
    suggestions = []
    total = sum(rejects.values())
    if total == 0:
        suggestions.append("No wire/OD combinations produced valid spring indices (C between 3 and 16). "
                           "Try widening the OD range.")
    else:
        if rejects["solid_too_tall"] > 0:
            suggestions.append(f"{rejects['solid_too_tall']} designs rejected: solid height exceeds compress-to. "
                               "Try increasing OD, reducing the solid height margin, or lowering target FPS/rate "
                               "(which reduces the required coil count).")
        if rejects["overstressed"] > 0:
            suggestions.append(f"{rejects['overstressed']} designs rejected: overstressed at compress-to. "
                               "Try increasing OD (allows thicker wire at safe stress), "
                               "increasing free length, or lowering target FPS/rate.")
        if rejects["spring_index"] > 0:
            suggestions.append(f"{rejects['spring_index']} designs rejected: spring index out of range (3-16). "
                               "Try widening the OD range.")
        if rejects["na_too_low"] > 0:
            suggestions.append(f"{rejects['na_too_low']} designs rejected: required coils < 1. "
                               "Try lowering target FPS/rate or increasing free length.")
    for s in suggestions:
        st.warning(s)
    st.stop()

# ══════════════════════════════════════════════
#  Results table
# ══════════════════════════════════════════════

n_total = len(candidates)
n_safe = sum(1 for c in candidates if c["safe_to_solid"])

_total_rejected = sum(rejects.values())
if _total_rejected > 0:
    _reject_labels = {
        "spring_index": "spring index out of range (3\u201316)",
        "solid_too_tall": "solid height exceeds compress-to",
        "overstressed": "overstressed at compress-to",
        "na_too_low": "required active coils < 1",
    }
    with st.expander(f"{_total_rejected} designs rejected"):
        for reason, count in rejects.items():
            if count > 0:
                st.write(f"- **{count}** {_reject_labels[reason]}")

if n_total > 10 and st.checkbox("Show only Pareto-optimal", value=True,
                                help="Keep only designs where no other design is better in both "
                                     "solid-height margin and stress utilization. Removes dominated options."):
    candidates = pareto_filter(candidates)
    st.caption(f"{len(candidates)} of {n_total} candidates (Pareto-optimal) — "
               f"{sum(1 for c in candidates if c['safe_to_solid'])} safe to solid")
else:
    st.caption(f"{n_total} candidates found — {n_safe} safe to solid")

df = pd.DataFrame([{
    "Wire d": format_wire_d(c["d"]),
    "OD": f'{c["OD"]:.3f}',
    "Na": f'{c["Na"]:g}',
    "Rate (lbf/in)": f'{c["k"]:.3f}',
    "Force @ CT (lbf)": f'{c["F_at_ct"]:.2f}',
    "Solid Ht (mm)": f'{c["Hs_mm"]:.1f}',
    "Margin (mm)": f'{c["margin_mm"]:.1f}',
    "Safe to Solid": "Yes" if c["safe_to_solid"] else "No",
    "Util @ CT": f'{c["util_at_ct"]:.0%}',
    "Util @ Solid": f'{c["util_at_solid"]:.0%}',
    "Index (C)": f'{c["C"]:.1f}',
} for c in candidates])

st.dataframe(
    df.style.apply(
        lambda row: [
            "background-color: rgba(34,197,94,0.12)" if row["Safe to Solid"] == "Yes"
            else "background-color: rgba(234,179,8,0.12)"
        ] * len(row),
        axis=1,
    ),
    use_container_width=True,
    hide_index=True,
    height=min(38 * len(df) + 40, 600),
)

st.download_button(
    "Export CSV", df.to_csv(index=False).encode("utf-8"),
    "spring_candidates.csv", "text/csv",
)

# ── Open in Analysis ──

st.markdown("**Open a design in Analysis**")
choice_idx = st.selectbox(
    "Select row", range(len(candidates)),
    format_func=lambda i: (
        f"#{i+1}: d={format_wire_d(candidates[i]['d'])}, OD={candidates[i]['OD']:.3f}, "
        f"Na={candidates[i]['Na']:g}, k={candidates[i]['k']:.3f} lbf/in"
        f" — {'safe to solid' if candidates[i]['safe_to_solid'] else 'safe at CT'}"
    ),
    label_visibility="collapsed",
)

sel = candidates[choice_idx]
from urllib.parse import urlencode
analysis_params = urlencode({
    "od": f"{sel['OD']:.3f}", "lf": f"{sel['Lf']:.2f}", "mat": wire_type,
    "d": f"{sel['d']:.6f}", "na": str(sel["Na"]), "end": end_type,
    "blaster": preset,
    "comp_from": f"{comp_from_mm:.1f}", "comp_to": f"{comp_to_mm:.1f}",
})
st.link_button("Open in Analysis", f"/?" + analysis_params)
