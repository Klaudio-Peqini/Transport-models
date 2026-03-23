from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Dict, Any

import json
import numpy as np
import pandas as pd


@dataclass(slots=True)
class PolymerBLParams:
    phi: float = 0.20
    L: float = 100.0
    area_m2: float = 1.0
    u: float = 1.0e-5
    Swi: float = 0.20
    Sor: float = 0.0
    nw: float = 2.0
    no: float = 2.0
    krw0: float = 0.80
    kro0: float = 1.00
    mu_w: float = 1.0
    mu_wp: float = 10.0
    mu_o: float = 5.0
    n_grid: int = 2001
    pvi_max: float = 2.0
    n_time: int = 500
    profile_pvi: tuple[float, ...] = (0.2, 0.5, 1.0)
    label: str = "polymer"

    @property
    def Sw_inj(self) -> float:
        return 1.0 - self.Sor


def effective_saturation(Sw: np.ndarray, Swi: float, Sor: float) -> np.ndarray:
    denom = max(1.0 - Swi - Sor, 1.0e-12)
    Se = (Sw - Swi) / denom
    return np.clip(Se, 0.0, 1.0)


def fractional_flow(Sw: np.ndarray, p: PolymerBLParams, water_viscosity: float) -> tuple[np.ndarray, np.ndarray]:
    Se = effective_saturation(Sw, p.Swi, p.Sor)
    krw = p.krw0 * Se**p.nw
    kro = p.kro0 * (1.0 - Se)**p.no
    lambda_w = krw / water_viscosity
    lambda_o = kro / p.mu_o
    fw = lambda_w / np.maximum(lambda_w + lambda_o, 1.0e-15)
    dfw = np.gradient(fw, Sw, edge_order=2)
    return fw, dfw


def shock_saturation(Sw: np.ndarray, fw: np.ndarray, dfw: np.ndarray, p: PolymerBLParams) -> float:
    fwi = np.interp(p.Swi, Sw, fw)
    mask = Sw > p.Swi + 1.0e-6
    secant = (fw[mask] - fwi) / (Sw[mask] - p.Swi)
    mismatch = np.abs(dfw[mask] - secant)
    idx_local = np.argmin(mismatch)
    return float(Sw[mask][idx_local])


def _case_data(p: PolymerBLParams, water_viscosity: float, label: str) -> Dict[str, Any]:
    Sw = np.linspace(p.Swi, p.Sw_inj, p.n_grid)
    fw, dfw = fractional_flow(Sw, p, water_viscosity)
    Swf = shock_saturation(Sw, fw, dfw, p)
    fwi = np.interp(p.Swi, Sw, fw)
    fwf = np.interp(Swf, Sw, fw)
    dfwf = np.interp(Swf, Sw, dfw)
    vD_shock = (fwf - fwi) / max(Swf - p.Swi, 1.0e-15)
    tD_bt = 1.0 / max(vD_shock, 1.0e-15)

    S_rare = np.linspace(Swf, p.Sw_inj, 1200)
    fw_rare = np.interp(S_rare, Sw, fw)
    dfdS_rare = np.interp(S_rare, Sw, dfw)
    order = np.argsort(dfdS_rare)
    vD_sorted = dfdS_rare[order]
    S_sorted = S_rare[order]
    # deduplicate for interpolation stability
    keep = np.concatenate(([True], np.diff(vD_sorted) > 1e-9))
    vD_sorted = vD_sorted[keep]
    S_sorted = S_sorted[keep]

    tD = np.linspace(1.0e-4, p.pvi_max, p.n_time)
    fw_prod = np.zeros_like(tD)
    Sw_prod = np.full_like(tD, p.Swi)
    for i, td in enumerate(tD):
        if td <= tD_bt:
            fw_prod[i] = fwi
            Sw_prod[i] = p.Swi
        else:
            target = 1.0 / td
            target = np.clip(target, float(vD_sorted.min()), float(vD_sorted.max()))
            Sw_out = np.interp(target, vD_sorted, S_sorted)
            Sw_prod[i] = Sw_out
            fw_prod[i] = np.interp(Sw_out, Sw, fw)

    dtD = np.gradient(tD)
    rf = np.cumsum((1.0 - fw_prod) * dtD) / max(1.0 - p.Swi - p.Sor, 1.0e-15)
    water_cut = fw_prod
    oil_cut = 1.0 - fw_prod

    profiles = []
    xD = np.linspace(0.0, 1.0, 600)
    for td in p.profile_pvi:
        prof = np.full_like(xD, p.Swi)
        if td <= 0:
            profiles.append(pd.DataFrame({"xD": xD, "Sw": prof, "PVI": td, "case": label}))
            continue
        xi = xD / td
        vmin = float(vD_sorted.min())
        vmax = float(vD_sorted.max())
        left = xi < vmin
        middle = (xi >= vmin) & (xi <= vD_shock)
        right = xi > vD_shock
        prof[left] = p.Sw_inj
        prof[middle] = np.interp(xi[middle], vD_sorted, S_sorted)
        prof[right] = p.Swi
        profiles.append(pd.DataFrame({"xD": xD, "Sw": prof, "PVI": td, "case": label}))
    profiles_df = pd.concat(profiles, ignore_index=True)

    area = p.area_m2
    pore_volume = p.phi * area * p.L
    ooip = pore_volume * (1.0 - p.Swi - p.Sor)
    q_total = p.u * area
    t_seconds = tD * p.phi * p.L / p.u
    cum_oil_m3 = rf * ooip
    water_rate_m3_s = water_cut * q_total
    oil_rate_m3_s = oil_cut * q_total

    curves = pd.DataFrame({
        "PVI": tD,
        "time_s": t_seconds,
        "water_cut": water_cut,
        "oil_cut": oil_cut,
        "Sw_prod": Sw_prod,
        "RF": rf,
        "cum_oil_m3": cum_oil_m3,
        "water_rate_m3_s": water_rate_m3_s,
        "oil_rate_m3_s": oil_rate_m3_s,
        "case": label,
    })

    fw_table = pd.DataFrame({"Sw": Sw, "fw": fw, "dfw_dSw": dfw, "case": label})
    summary = {
        "case": label,
        "mu_water_cP": water_viscosity,
        "Swi": p.Swi,
        "Sw_inj": p.Sw_inj,
        "Sw_shock": Swf,
        "fw_shock": fwf,
        "vD_shock": vD_shock,
        "PVI_breakthrough": tD_bt,
        "RF_at_breakthrough": float(np.interp(tD_bt, tD, rf)),
        "RF_final": float(rf[-1]),
        "water_cut_final": float(water_cut[-1]),
    }
    return {"params": asdict(p), "summary": summary, "fractional_flow": fw_table, "curves": curves, "profiles": profiles_df}


def compute_case(params: PolymerBLParams, use_polymer: bool = True) -> Dict[str, Any]:
    mu = params.mu_wp if use_polymer else params.mu_w
    label = "polymer" if use_polymer else "waterflood"
    return _case_data(params, mu, label)


def compare_cases(params: PolymerBLParams) -> Dict[str, Any]:
    water = compute_case(params, use_polymer=False)
    polymer = compute_case(params, use_polymer=True)
    return {
        "waterflood": water,
        "polymer": polymer,
        "summaries": pd.DataFrame([water["summary"], polymer["summary"]]),
        "fractional_flow": pd.concat([water["fractional_flow"], polymer["fractional_flow"]], ignore_index=True),
        "curves": pd.concat([water["curves"], polymer["curves"]], ignore_index=True),
        "profiles": pd.concat([water["profiles"], polymer["profiles"]], ignore_index=True),
    }


def run_transport_bridge(params: PolymerBLParams, wwtp_results: pd.DataFrame | dict) -> pd.DataFrame:
    comparison = compare_cases(params)
    curves = comparison["curves"].copy()
    polymer_prod = curves[curves["case"] == "polymer"][["PVI", "water_cut"]].rename(columns={"water_cut": "normalized_breakthrough"})
    water_prod = curves[curves["case"] == "waterflood"][["PVI", "water_cut"]].rename(columns={"water_cut": "normalized_breakthrough"})
    polymer_prod["domain"] = "Polymer front"
    water_prod["domain"] = "Waterflood front"

    if isinstance(wwtp_results, dict):
        if "timeseries" not in wwtp_results:
            raise ValueError("When a dict is provided, it must contain a 'timeseries' DataFrame")
        ts = wwtp_results["timeseries"].copy()
    else:
        ts = wwtp_results.copy()

    time_col = "t_h" if "t_h" in ts.columns else ("time_h" if "time_h" in ts.columns else None)
    if time_col is None:
        raise ValueError("wwtp_results must contain either 't_h' or 'time_h'")
    max_t = float(ts[time_col].max())
    bridge_rows = []
    for species in ["CBZ_OUT", "DCF_OUT", "CBZ_gac", "DCF_gac"]:
        if species in ts.columns:
            pretty = species.replace("_OUT", "").replace("_gac", "")
            df = pd.DataFrame({
                "PVI": ts[time_col] / max_t,
                "normalized_breakthrough": ts[species] / max(ts[species].max(), 1e-15),
                "domain": f"Wastewater {pretty} outlet",
            })
            bridge_rows.append(df)
    bridge = pd.concat([polymer_prod, water_prod] + bridge_rows, ignore_index=True)
    return bridge


def save_case_bundle(bundle: Dict[str, Any], output_dir: str | Path) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    bundle["summaries"].to_csv(out / "polymer_summary.csv", index=False)
    bundle["fractional_flow"].to_csv(out / "fractional_flow_table.csv", index=False)
    bundle["curves"].to_csv(out / "polymer_production_curves.csv", index=False)
    bundle["profiles"].to_csv(out / "polymer_profiles.csv", index=False)
    with open(out / "polymer_config.json", "w", encoding="utf-8") as f:
        json.dump(bundle["polymer"]["params"], f, indent=2)
