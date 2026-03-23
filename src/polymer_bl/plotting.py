from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _ensure_dir(output_dir: str | Path) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out


def plot_fractional_flow(bundle: Dict[str, Any], output_dir: str | Path) -> None:
    out = _ensure_dir(output_dir)
    ff = bundle["fractional_flow"]
    summaries = bundle["summaries"].set_index("case")
    plt.figure(figsize=(8.2, 5.4))
    for case, grp in ff.groupby("case"):
        plt.plot(grp["Sw"], grp["fw"], linewidth=2.2, label=f"{case}: $f_w(S_w)$")
        Swf = summaries.loc[case, "Sw_shock"]
        fwf = summaries.loc[case, "fw_shock"]
        slope = summaries.loc[case, "vD_shock"]
        xline = np.array([summaries.loc[case, "Swi"], Swf])
        yline = slope * (xline - summaries.loc[case, "Swi"])
        plt.plot(xline, yline, linestyle="--", linewidth=1.6, label=f"Tangjenta Welge ({case})")
        plt.axvline(Swf, linestyle=":", linewidth=1.3)
        plt.scatter([Swf], [fwf], s=32)
    plt.xlabel(r"Saturimi i ujit $S_w$")
    plt.ylabel(r"Rrjedha fraksionare e ujit $f_w$")
    plt.title("Kurba e rrjedhës fraksionare: waterflood kundrejt polymer flood")
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True, fontsize=9)
    plt.tight_layout()
    plt.savefig(out / "fractional_flow_comparison.png", dpi=300)
    plt.close()


def plot_profiles(bundle: Dict[str, Any], output_dir: str | Path) -> None:
    out = _ensure_dir(output_dir)
    profiles = bundle["profiles"]
    pvls = sorted(profiles["PVI"].unique())
    fig, axes = plt.subplots(1, len(pvls), figsize=(5.2 * len(pvls), 4.8), sharey=True)
    if len(pvls) == 1:
        axes = [axes]
    for ax, pvi in zip(axes, pvls):
        sub = profiles[profiles["PVI"] == pvi]
        for case, grp in sub.groupby("case"):
            ax.plot(grp["xD"], grp["Sw"], linewidth=2.2, label=case)
        ax.set_title(f"Profili i saturimit në PVI = {pvi:.2f}")
        ax.set_xlabel(r"Pozicioni i reduktuar $x_D=x/L$")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"Saturimi i ujit $S_w$")
    axes[-1].legend(frameon=True)
    plt.tight_layout()
    plt.savefig(out / "saturation_profiles_comparison.png", dpi=300)
    plt.close(fig)


def plot_production(bundle: Dict[str, Any], output_dir: str | Path) -> None:
    out = _ensure_dir(output_dir)
    curves = bundle["curves"]
    summaries = bundle["summaries"].set_index("case")
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    for case, grp in curves.groupby("case"):
        axes[0].plot(grp["PVI"], grp["water_cut"], linewidth=2.2, label=case)
        axes[1].plot(grp["PVI"], grp["RF"], linewidth=2.2, label=case)
        axes[0].axvline(summaries.loc[case, "PVI_breakthrough"], linestyle="--", linewidth=1.2)
    axes[0].set_xlabel("Pore volumes injected, PVI [-]")
    axes[0].set_ylabel("Water cut në prodhim [-]")
    axes[0].set_title("Kurba e breakthrough në prodhues")
    axes[0].grid(True, alpha=0.3)
    axes[1].set_xlabel("Pore volumes injected, PVI [-]")
    axes[1].set_ylabel("Faktori i rikuperimit RF [-]")
    axes[1].set_title("Rikuperimi kumulativ i naftës")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(frameon=True)
    plt.tight_layout()
    plt.savefig(out / "production_curves_comparison.png", dpi=300)
    plt.close(fig)


def plot_transport_bridge(bridge_df: pd.DataFrame, output_dir: str | Path) -> None:
    out = _ensure_dir(output_dir)
    plt.figure(figsize=(8.8, 5.2))
    for domain, grp in bridge_df.groupby("domain"):
        plt.plot(grp["PVI"], grp["normalized_breakthrough"], linewidth=2.0, label=domain)
    plt.xlabel("Koha / vëllimi i injektuar i reduktuar [-]")
    plt.ylabel("Breakthrough i normalizuar [-]")
    plt.title("Lidhja konceptuale midis frontit Buckley–Leverett dhe breakthrough në modelin e ujërave")
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True, fontsize=9)
    plt.tight_layout()
    plt.savefig(out / "transport_bridge.png", dpi=300)
    plt.close()
