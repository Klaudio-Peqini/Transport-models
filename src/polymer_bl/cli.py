from __future__ import annotations

import argparse
from pathlib import Path

from .model import PolymerBLParams, compare_cases, save_case_bundle, run_transport_bridge
from .plotting import plot_fractional_flow, plot_profiles, plot_production, plot_transport_bridge


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Analytical Buckley-Leverett polymer flooding toolkit")
    p.add_argument("--phi", type=float, default=0.20)
    p.add_argument("--L", type=float, default=100.0)
    p.add_argument("--u", type=float, default=1.0e-5)
    p.add_argument("--Swi", type=float, default=0.20)
    p.add_argument("--Sor", type=float, default=0.0)
    p.add_argument("--nw", type=float, default=2.0)
    p.add_argument("--no", type=float, default=2.0)
    p.add_argument("--krw0", type=float, default=0.80)
    p.add_argument("--kro0", type=float, default=1.00)
    p.add_argument("--mu-w", dest="mu_w", type=float, default=1.0)
    p.add_argument("--mu-wp", dest="mu_wp", type=float, default=10.0)
    p.add_argument("--mu-o", dest="mu_o", type=float, default=5.0)
    p.add_argument("--pvi-max", type=float, default=2.0)
    p.add_argument("--n-grid", type=int, default=2001)
    p.add_argument("--n-time", type=int, default=500)
    p.add_argument("--profile-pvi", type=float, nargs="*", default=[0.2, 0.5, 1.0])
    p.add_argument("--output-dir", default="outputs_polymer")
    return p


def main() -> None:
    args = build_parser().parse_args()
    params = PolymerBLParams(
        phi=args.phi, L=args.L, u=args.u, Swi=args.Swi, Sor=args.Sor,
        nw=args.nw, no=args.no, krw0=args.krw0, kro0=args.kro0,
        mu_w=args.mu_w, mu_wp=args.mu_wp, mu_o=args.mu_o,
        pvi_max=args.pvi_max, n_grid=args.n_grid, n_time=args.n_time,
        profile_pvi=tuple(args.profile_pvi),
    )
    bundle = compare_cases(params)
    out = Path(args.output_dir)
    save_case_bundle(bundle, out)
    plot_fractional_flow(bundle, out)
    plot_profiles(bundle, out)
    plot_production(bundle, out)
    print(f"Saved polymer Buckley-Leverett outputs to: {out}")


if __name__ == "__main__":
    main()
