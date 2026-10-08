"""
LES Backward-Facing Step — RANS Probes Post-processing Pipeline
===============================================================
Horses3d high-order LES · 700 snapshots/CTU
Probe layout : 14x × 10y × 100z = 14 000 probes per snapshot
Variables    : rho, u, v, w, p   (names as stored in the HDF5 file)
U_REF = 1,  H = 1  (fixed)

Input
-----
A single HDF5 file written by Horses3d Monitor_WriteFileProbesHDF5:

  /coordinates  (3, nProbes_total)  float64  — x/y/z for ALL probes
  /time         (n_steps,)          float64
  /iteration    (n_steps,)          int32
  /<varname>    (n_steps, nProbes_total) float64   one dataset per variable

The first N_RANS_PROBES rows correspond to the RANS probe layout
(14x × 10y × 100z = 14 000).

Key geometry note
-----------------
The 10 wall-normal (y) probe positions are LOCAL to each x station:
each x has its own y grid that follows the step geometry.
y_locs is therefore a 2-D array of shape (nx, ny).
z is shared and uniform across all stations.

Analyses
--------
  1. Mean & rms profiles         one PNG per x-section (2 rows × 4 cols)
  2. Streamwise evolution        one PNG, 8 stacked panels (waterfall)
  3. Spanwise energy spectra     one PNG per x-section (all y)
  4. Two-point correlations      one PNG with all selected stations
  5. Temporal power spectra      one PNG per selected station

Output structure
----------------
  POSTPRO/Probes_RANS/
    data/
      mean_profiles.npz
      spectra_spanwise.npz
      two_point_corr.npz
      spectra_temporal.npz
    figures/
      profiles_x??.png
      profiles_streamwise_evolution.png
      spectra_spanwise_x??.png
      tpc_all_stations.png
      spectra_temporal_x??_y??.png

Usage
-----
  python plot_probes_RANS_postpro.py --h5-file RESULTS/run.probes.h5 \\
                                     [--n-rans 14000] \\
                                     [--nx 14] [--ny 10] [--nz 100] \\
                                     [--var-u u] [--var-v v] [--var-w w] \\
                                     [--var-p pressure] [--var-rho rho] \\
                                     [--output-dir POSTPRO] [--dt-ctu 0.001429]
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import welch

# ── Fixed physical parameters ─────────────────────────────────────────────────
U_REF = 1.0
H     = 1.0

# TPC & temporal spectra stations: even ix indices (0,2,4,…)
# IY_TARGET: wall-normal index; None → ny//2
IY_TARGET = None

# Streamwise evolution plot x-axis range
X_MIN, X_MAX = -1.0, 3.0

# Number of RANS probes (first N rows in the HDF5 file)
N_RANS_PROBES = 14_000


# ── HDF5 reader ───────────────────────────────────────────────────────────────

def load_h5(h5_file: Path, nx: int, ny: int, nz: int,
            n_rans: int, var_names: dict):
    """
    Read the Horses3d .probes.h5 file and return the RANS sub-set.

    Parameters
    ----------
    h5_file   : path to the .probes.h5 file
    nx,ny,nz  : RANS grid dimensions (nx*ny*nz must equal n_rans)
    n_rans    : number of RANS probes (first rows in the file)
    var_names : mapping from logical name to HDF5 dataset name
                keys: "u","v","w","p","rho"

    Returns
    -------
    x_locs     : (nx,)       x coordinates
    y_locs     : (nx, ny)    y coordinates (wall-adapted per x station)
    z_locs     : (nz,)       z coordinates
    probe_order: list of (ix, iy, iz) for each of the n_rans probe rows
    G          : (n_steps, nx, ny, nz, 4)  u v w p
    time       : (n_steps,)  simulation time
    iteration  : (n_steps,)  iteration number
    """
    if nx * ny * nz != n_rans:
        sys.exit(f"[ERROR] nx({nx}) × ny({ny}) × nz({nz}) = {nx*ny*nz} ≠ n_rans({n_rans})")

    with h5py.File(h5_file, "r") as fh:
        # ── coordinates (RANS subset only) ────────────────────────────────────
        # shape in file: (3, nProbes_total) — Fortran column-major
        coords = fh["coordinates"][:, :n_rans]   # (3, n_rans)
        x_all, y_all, z_all = coords[0], coords[1], coords[2]

        x_locs = np.unique(x_all)
        z_locs = np.unique(z_all)
        if len(x_locs) != nx:
            sys.exit(f"[ERROR] Found {len(x_locs)} unique x values, expected {nx}")
        if len(z_locs) != nz:
            sys.exit(f"[ERROR] Found {len(z_locs)} unique z values, expected {nz}")

        y_per_x = [np.unique(y_all[x_all == xv]) for xv in x_locs]
        ny_vals = [len(y) for y in y_per_x]
        if len(set(ny_vals)) != 1 or ny_vals[0] != ny:
            sys.exit(f"[ERROR] Unexpected ny distribution across x: {ny_vals}")
        y_locs = np.array(y_per_x)   # (nx, ny)

        probe_order = []
        for xv, yv, zv in zip(x_all, y_all, z_all):
            ix = int(np.searchsorted(x_locs, xv))
            iy = int(np.searchsorted(y_locs[ix], yv))
            iz = int(np.searchsorted(z_locs, zv))
            probe_order.append((ix, iy, iz))

        # ── time / iteration ──────────────────────────────────────────────────
        time      = fh["time"][:]
        iteration = fh["iteration"][:]
        n_steps   = len(time)
        print(f"  Steps    : {n_steps}")
        print(f"  Time     : {time[0]:.6g} → {time[-1]:.6g}")
        print(f"  Iteration: {iteration[0]} → {iteration[-1]}")

        # ── flow variables  ───────────────────────────────────────────────────
        # Each dataset: (n_steps, nProbes_total) — read RANS columns only
        G = np.empty((n_steps, nx, ny, nz, 4), dtype=np.float64)
        for vi, key in enumerate(["u", "v", "w", "p"]):
            dset_name = var_names[key]
            if dset_name not in fh:
                sys.exit(f"[ERROR] Dataset '{dset_name}' not found in {h5_file}. "
                         f"Available: {list(fh.keys())}")
            raw = fh[dset_name][:, :n_rans]   # (n_steps, n_rans)
            print(f"  Loading  : {dset_name}  {raw.shape}")
            for p_idx, (ix, iy, iz) in enumerate(probe_order):
                G[:, ix, iy, iz, vi] = raw[:, p_idx]

    return x_locs, y_locs, z_locs, probe_order, G, time, iteration


# ── Analysis 1 — Statistics ───────────────────────────────────────────────────

def compute_statistics(G: np.ndarray) -> dict:
    """Time-and-spanwise statistics, shape (nx, ny)."""
    mean_t  = G.mean(axis=0)
    mean_xy = mean_t.mean(axis=2)
    fluct   = G - mean_t[np.newaxis]
    rms_xy  = np.sqrt((fluct ** 2).mean(axis=(0, 3)))
    uv      = (fluct[..., 0] * fluct[..., 1]).mean(axis=(0, 3))
    return dict(
        mean_u=mean_xy[..., 0], mean_v=mean_xy[..., 1],
        mean_w=mean_xy[..., 2], mean_p=mean_xy[..., 3],
        rms_u =rms_xy[..., 0],  rms_v =rms_xy[..., 1],
        rms_w =rms_xy[..., 2],  rms_p =rms_xy[..., 3],
        uv_stress=-uv,
    )


# ── Analysis 1 — Profiles per x-section ──────────────────────────────────────

def plot_profiles(stats: dict, x_locs, y_locs, out_dir: Path):
    """2 rows × 4 cols per x: mean (top) + rms (bottom)."""
    fig_dir    = out_dir / "figures"
    mean_items = [("mean_u","⟨u⟩"),("mean_v","⟨v⟩"),("mean_w","⟨w⟩"),("mean_p","⟨p⟩")]
    rms_items  = [("rms_u","u_rms"),("rms_v","v_rms"),("rms_w","w_rms"),("uv_stress","−⟨u′v′⟩")]

    for ix, x_val in enumerate(x_locs):
        y_n = y_locs[ix] / H
        fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharey=True)
        fig.suptitle(f"Profiles — x/H = {x_val/H:.2f}", fontsize=12)

        for col, (key, label) in enumerate(mean_items):
            ax = axes[0, col]
            ax.plot(stats[key][ix], y_n, "o-", ms=4, lw=1.2, color="tab:blue")
            ax.axvline(0, color="gray", lw=0.7, ls="--")
            ax.set_xlabel(label, fontsize=10)
            ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
            ax.grid(True, lw=0.4, alpha=0.5)
            if col == 0: ax.set_ylabel("y", fontsize=10)

        for col, (key, label) in enumerate(rms_items):
            ax = axes[1, col]
            ax.plot(stats[key][ix], y_n, "s-", ms=4, lw=1.2, color="tab:orange")
            ax.axvline(0, color="gray", lw=0.7, ls="--")
            ax.set_xlabel(label, fontsize=10)
            ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
            ax.grid(True, lw=0.4, alpha=0.5)
            if col == 0: ax.set_ylabel("y", fontsize=10)

        fig.tight_layout()
        fig.savefig(fig_dir / f"profiles_x{ix:02d}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"  [profiles] → {fig_dir}")


# ── Analysis 2 — Streamwise evolution (waterfall) ─────────────────────────────

def _draw_geometry(axes, y_locs):
    """Draw rounded-step lower wall on every axis."""
    x_geo  = np.linspace(X_MIN, X_MAX, 4000)
    y_wall = np.where(x_geo <= 0, 1.0,
             np.where(x_geo <= 1.0,
                      1.0 - 10.0*x_geo**3 + 15.0*x_geo**4 - 6.0*x_geo**5,
                      0.0))

    ix_up = 0
    ix_dn = -1
    H_phys       = float(y_locs[ix_up].min())
    y_floor_phys = float(y_locs[ix_dn].min())
    y_wall_phys  = y_floor_phys + y_wall * (H_phys - y_floor_phys)
    y_fill_floor = y_floor_phys - 0.005

    for ax in axes:
        ax.plot(x_geo, y_wall_phys, color="black", lw=1.8, zorder=6)
        ax.fill_between(x_geo, y_wall_phys, y_fill_floor,
                        color="lightgray", alpha=0.6, zorder=5)
        ax.set_ylim(y_fill_floor, y_locs.max() + 0.01 * (y_locs.max() - y_locs.min()))
        ax.set_xlabel("x/H", fontsize=9)
        ax.tick_params(labelbottom=True, labelsize=8)


def plot_streamwise_evolution(stats: dict, x_locs, y_locs, out_dir: Path):
    """8 stacked panels (one per variable), waterfall offset profiles."""
    fig_dir = out_dir / "figures"
    variables = [
        ("mean_u","⟨u⟩"), ("mean_v","⟨v⟩"), ("mean_w","⟨w⟩"), ("mean_p","⟨p⟩"),
        ("rms_u","u_rms"), ("rms_v","v_rms"), ("rms_w","w_rms"), ("uv_stress","−⟨u′v′⟩"),
    ]

    x_gaps        = np.diff(x_locs / H)
    PROFILE_WIDTH = float(np.min(x_gaps)) * 0.8
    cmap          = plt.get_cmap("viridis")
    colors        = [cmap(i / max(len(x_locs)-1, 1)) for i in range(len(x_locs))]

    y_min  = y_locs.min()
    y_max  = y_locs.max()
    y_span = y_max - y_min
    x_span = X_MAX - X_MIN
    FIG_W  = 14.0
    FIG_H  = max(FIG_W * (y_span / x_span), 1.5)

    n_vars   = len(variables)
    GAP      = 0.05
    canvas_h = FIG_H * n_vars + GAP * (n_vars - 1) + 0.6
    fig      = plt.figure(figsize=(FIG_W, canvas_h))

    left, right  = 0.07, 0.97
    top_pad      = 0.6 / canvas_h
    bot_pad      = 0.03
    ax_h_frac    = FIG_H / canvas_h
    gap_frac     = GAP  / canvas_h

    axes = []
    for i in range(n_vars):
        bottom = 1.0 - top_pad - (i + 1) * ax_h_frac - i * gap_frac
        axes.append(fig.add_axes([left, bottom, right - left, ax_h_frac]))

    for ax, (key, label) in zip(axes, variables):
        field   = stats[key]
        centred = np.array([field[ix] - field[ix].mean() for ix in range(len(x_locs))])
        scale   = np.ptp(centred)
        scale   = scale if scale > 1e-30 else 1.0

        for ix, x_val in enumerate(x_locs):
            x_base = x_val / H
            if x_base < X_MIN or x_base > X_MAX:
                continue
            profile = (centred[ix] / scale) * PROFILE_WIDTH
            ax.plot(x_base + profile, y_locs[ix],
                    lw=1.3, color=colors[ix], label=f"x/H={x_val/H:.2f}")
            if key != "mean_p":
                ax.axvline(x_base, color=colors[ix], lw=0.5, ls="--", alpha=0.35)

        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(y_min - 0.01*y_span, y_max + 0.01*y_span)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("y", fontsize=9)
        var_short = label.replace("⟨","").replace("⟩","").replace("−","-").replace("′","'")
        ax.set_title(f"{label}   [scale = {var_short}_max - {var_short}_min = {scale:.3g}]",
                     fontsize=9, loc="left")
        ax.grid(True, lw=0.3, alpha=0.4)
        ax.tick_params(labelsize=8)

    _draw_geometry(axes, y_locs)

    handles, labels_leg = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_leg, loc="upper center", ncol=len(x_locs),
               fontsize=7, frameon=True, bbox_to_anchor=(0.5, 1.0 - top_pad / 2))

    fname = fig_dir / "profiles_streamwise_evolution.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [evolution] → {fname}")


# ── Analysis 3 — Spanwise energy spectra ─────────────────────────────────────

def compute_spanwise_spectra(G: np.ndarray, z_locs: np.ndarray) -> dict:
    dz = z_locs[1] - z_locs[0]
    nz = len(z_locs)
    kz = np.fft.rfftfreq(nz, d=dz)
    E  = {}
    for vi, var in enumerate(["u","v","w"]):
        fft_z = np.fft.rfft(G[..., vi], axis=3)
        psd   = (np.abs(fft_z)**2) / nz * 2
        psd[:, :, :, 0] /= 2
        if nz % 2 == 0: psd[:, :, :, -1] /= 2
        E[f"E_{var}{var}"] = psd.mean(axis=0)
    return dict(kz=kz, **E)


def plot_spanwise_spectra(spec: dict, x_locs, y_locs, out_dir: Path):
    fig_dir = out_dir / "figures"
    kz, mask = spec["kz"], spec["kz"] > 0
    ny    = y_locs.shape[1]
    ncols = (ny + 1) // 2

    for ix, x_val in enumerate(x_locs):
        fig, axes = plt.subplots(2, ncols, figsize=(4*ncols, 8), sharey=True)
        axes = np.array(axes).reshape(2, ncols)
        fig.suptitle(f"Spanwise spectra — x/H = {x_val/H:.2f}", fontsize=11)

        for iy in range(ny):
            row, col = divmod(iy, ncols)
            ax = axes[row, col]
            for var, color, ls in [("uu","tab:blue","-"),("vv","tab:orange","--"),("ww","tab:green","-.")]:
                ax.loglog(kz[mask]*H, spec[f"E_{var}"][ix,iy,mask], lw=1.2, color=color, ls=ls, label=f"$E_{{{var}}}$")
            k_ref = kz[mask][len(kz[mask])//3]
            E_ref = spec["E_uu"][ix, iy, len(kz[mask])//3]
            ax.loglog(kz[mask]*H, E_ref*(kz[mask]/k_ref)**(-5/3), "k--", lw=0.8, label="$-5/3$")
            ax.set_title(f"y={y_locs[ix,iy]:.4f}", fontsize=9)
            ax.set_xlabel("$k_z H$", fontsize=9)
            if col == 0: ax.set_ylabel("$E(k_z)$", fontsize=9)
            ax.legend(fontsize=7)
            ax.grid(True, which="both", lw=0.3, alpha=0.4)

        for iy in range(ny, 2*ncols):
            row, col = divmod(iy, ncols)
            axes[row, col].set_visible(False)

        fig.tight_layout()
        fig.savefig(fig_dir / f"spectra_spanwise_x{ix:02d}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"  [spanwise spectra] → {fig_dir}")


# ── Analysis 4 — Two-point correlations ──────────────────────────────────────

def resolve_stations(nx: int, ny: int):
    iy = ny // 2 if IY_TARGET is None else int(IY_TARGET)
    return [(ix, iy) for ix in range(0, nx, 2)]


def compute_two_point_correlation(G: np.ndarray, stations: list) -> dict:
    nz, results = G.shape[3], {}
    for (ix, iy) in stations:
        sub = {"dz_idx": np.arange(nz)}
        for vi, var in enumerate(["u","v","w","p"]):
            signal = G[:, ix, iy, :, vi]
            fluct  = signal - signal.mean(axis=0)
            fluct  = fluct  - fluct.mean(axis=1, keepdims=True)
            var0   = (fluct**2).mean()
            if var0 < 1e-30:
                sub[f"R_{var}{var}"] = np.ones(nz)
                continue
            F    = np.fft.rfft(fluct, axis=1)
            corr = np.fft.irfft(F * np.conj(F), n=nz, axis=1)
            sub[f"R_{var}{var}"] = corr.mean(axis=0) / (nz * var0)
        results[(ix, iy)] = sub
    return results


def plot_two_point_correlations(tpc: dict, z_locs, x_locs, y_locs, out_dir: Path):
    fig_dir  = out_dir / "figures"
    dz, nz   = z_locs[1] - z_locs[0], len(z_locs)
    half     = nz // 2
    stations = list(tpc.keys())
    n_stat   = len(stations)

    fig, axes = plt.subplots(n_stat, 4, figsize=(16, 3.5*n_stat), squeeze=False)
    fig.suptitle("Two-point correlations $R(\\Delta z)$", fontsize=12)

    for row, (ix, iy) in enumerate(stations):
        lag_z = np.arange(half) * dz / H
        for col, (var, color) in enumerate([("uu","tab:blue"),("vv","tab:orange"),
                                             ("ww","tab:green"),("pp","tab:red")]):
            ax = axes[row, col]
            R  = tpc[(ix,iy)][f"R_{var}"][:half]
            ax.plot(lag_z, R, lw=1.4, color=color)
            ax.axhline(0, color="gray", lw=0.7, ls="--")
            if var == "uu":
                zeros = np.where(np.diff(np.sign(R)))[0]
                if len(zeros):
                    iz0  = zeros[0]
                    L_uu = float(np.trapz(R[:iz0+1], lag_z[:iz0+1]))
                    ax.axvline(lag_z[iz0], color=color, lw=0.8, ls=":",
                               label=f"$L_{{uu}}\\approx{L_uu:.3f}H$")
                    ax.legend(fontsize=8)
            ax.set_xlim(0, lag_z[-1])
            ax.set_xlabel("$\\Delta z/H$", fontsize=9)
            ax.set_ylabel(f"$R_{{{var}}}$", fontsize=9)
            ax.grid(True, lw=0.4, alpha=0.5)
            if col == 0:
                ax.set_title(f"x/H={x_locs[ix]/H:.2f}  y={y_locs[ix,iy]:.4f}",
                             fontsize=9, loc="left")

    fig.tight_layout()
    fig.savefig(fig_dir / "tpc_all_stations.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [two-point] → {fig_dir / 'tpc_all_stations.png'}")


# ── Analysis 5 — Temporal power spectra ──────────────────────────────────────

def compute_temporal_spectra(G: np.ndarray, stations: list, dt_ctu: float) -> dict:
    n_snap  = G.shape[0]
    nperseg = max(64, n_snap // 4)
    fs      = 1.0 / dt_ctu
    results = {}
    for (ix, iy) in stations:
        sub = {}
        for vi, var in enumerate(["u","v","w","p"]):
            signal = G[:, ix, iy, :, vi].mean(axis=1)
            signal = signal - signal.mean()
            freq, psd = welch(signal, fs=fs, nperseg=nperseg, window="hann", scaling="density")
            sub[f"PSD_{var}{var}"] = psd
        sub["freq"] = freq
        results[(ix, iy)] = sub
    return results


def plot_temporal_spectra(tspec: dict, x_locs, y_locs, out_dir: Path):
    fig_dir = out_dir / "figures"
    for (ix, iy), sub in tspec.items():
        freq, mask = sub["freq"], sub["freq"] > 0
        fig, ax = plt.subplots(figsize=(7, 4))
        for var, color, ls in [("uu","tab:blue","-"),("vv","tab:orange","--"),
                                ("ww","tab:green","-."),("pp","tab:red",":")]:
            ax.loglog(freq[mask], sub[f"PSD_{var}"][mask], lw=1.2, color=color, ls=ls,
                      label=f"$PSD_{{{var}}}$")
        imid  = len(freq[mask]) // 3
        f_ref = freq[mask][imid]
        p_ref = sub["PSD_uu"][mask][imid]
        ax.loglog(freq[mask], p_ref*(freq[mask]/f_ref)**(-5/3), "k--", lw=0.8, label="$-5/3$")
        ax.set_xlabel("$f$ [CTU$^{-1}$]", fontsize=10)
        ax.set_ylabel("PSD", fontsize=10)
        ax.set_title(f"Temporal spectrum  x/H={x_locs[ix]/H:.2f}  y={y_locs[ix,iy]:.4f}",
                     fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, which="both", lw=0.3, alpha=0.4)
        fig.tight_layout()
        fig.savefig(fig_dir / f"spectra_temporal_x{ix:02d}_y{iy:02d}.png",
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(f"  [temporal spectra] → {fig_dir}")


# ── Save ──────────────────────────────────────────────────────────────────────

def save_results(out_dir, stats, spec_z, tpc, tspec, x_locs, y_locs, z_locs):
    data_dir = out_dir / "data"
    coords   = dict(x_locs=x_locs, y_locs=y_locs, z_locs=z_locs)
    np.savez_compressed(data_dir / "mean_profiles.npz",    **coords, **stats)
    np.savez_compressed(data_dir / "spectra_spanwise.npz", **coords, **spec_z)

    def flatten(d):
        out = {}
        for (ix, iy), sub in d.items():
            pfx = f"x{ix:02d}_y{iy:02d}"
            for k, v in sub.items():
                out[f"{pfx}_{k}"] = np.asarray(v)
        return out

    np.savez_compressed(data_dir / "two_point_corr.npz",   **coords, **flatten(tpc))
    np.savez_compressed(data_dir / "spectra_temporal.npz", **coords, **flatten(tspec))
    print(f"  [data] → {data_dir}")


# ── Main ──────────────────────────────────────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(description="RANS probes post-processing pipeline (HDF5 input)")
    p.add_argument("--h5-file",    required=True,
                   help="Horses3d HDF5 probe file, e.g. run.probes.h5")
    p.add_argument("--n-rans",     type=int, default=N_RANS_PROBES,
                   help="Number of RANS probes (first N rows in the HDF5). Default: 14000")
    p.add_argument("--nx",         type=int, default=14,
                   help="Number of x stations. Default: 14")
    p.add_argument("--ny",         type=int, default=10,
                   help="Number of wall-normal points per station. Default: 10")
    p.add_argument("--nz",         type=int, default=100,
                   help="Number of spanwise points. Default: 100")
    p.add_argument("--var-u",      default="u",
                   help="HDF5 dataset name for u velocity. Default: u")
    p.add_argument("--var-v",      default="v",
                   help="HDF5 dataset name for v velocity. Default: v")
    p.add_argument("--var-w",      default="w",
                   help="HDF5 dataset name for w velocity. Default: w")
    p.add_argument("--var-p",      default="pressure",
                   help="HDF5 dataset name for pressure. Default: pressure")
    p.add_argument("--var-rho",    default="rho",
                   help="HDF5 dataset name for density (unused in plots, kept for save). Default: rho")
    p.add_argument("--output-dir", default="POSTPRO",
                   help="Root output dir; results go into <output-dir>/Probes_RANS/")
    p.add_argument("--dt-ctu",     type=float, default=1.0/700.0,
                   help="Time between snapshots in CTU. Default: 1/700")
    p.add_argument("--skip-plots", action="store_true")
    return p


def run(args):
    h5_file = Path(args.h5_file)
    if not h5_file.exists():
        sys.exit(f"[ERROR] File not found: {h5_file}")

    out_dir = Path(args.output_dir) / "Probes_RANS"
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    (out_dir / "data").mkdir(parents=True, exist_ok=True)

    var_names = dict(u=args.var_u, v=args.var_v, w=args.var_w,
                     p=args.var_p, rho=args.var_rho)

    # ── 0. Load HDF5 ──────────────────────────────────────────────────────
    print(f"\n[0/5] Reading HDF5: {h5_file} …")
    print(f"  RANS probes: first {args.n_rans} of all probes")
    x_locs, y_locs, z_locs, probe_order, G, time, iteration = load_h5(
        h5_file, args.nx, args.ny, args.nz, args.n_rans, var_names)

    n_steps = G.shape[0]
    nx, ny, nz = args.nx, args.ny, args.nz
    print(f"  Grid : {nx}x × {ny}y × {nz}z = {nx*ny*nz} RANS probes")
    print(f"  x/H  : {[f'{v/H:.2f}' for v in x_locs]}")
    print(f"  Lz   : {z_locs[0]:.4f} → {z_locs[-1]:.4f}   Δz={z_locs[1]-z_locs[0]:.5f}")
    print(f"  G    : {G.shape}  ({G.nbytes/1e9:.2f} GB)")

    stations = resolve_stations(nx, ny)
    print(f"  Stations (even ix, iy={stations[0][1]}):")
    for ix, iy in stations:
        print(f"    ix={ix}  x/H={x_locs[ix]/H:.2f}  y={y_locs[ix,iy]:.4f}")

    # ── 1. Statistics & profiles ──────────────────────────────────────────
    print("\n[1/5] Statistics …")
    stats = compute_statistics(G)
    if not args.skip_plots:
        plot_profiles(stats, x_locs, y_locs, out_dir)
        plot_streamwise_evolution(stats, x_locs, y_locs, out_dir)

    # ── 2. Spanwise spectra ───────────────────────────────────────────────
    print("\n[2/5] Spanwise spectra …")
    spec_z = compute_spanwise_spectra(G, z_locs)
    if not args.skip_plots:
        plot_spanwise_spectra(spec_z, x_locs, y_locs, out_dir)

    # ── 3. Two-point correlations ─────────────────────────────────────────
    print("\n[3/5] Two-point correlations …")
    tpc = compute_two_point_correlation(G, stations)
    if not args.skip_plots:
        plot_two_point_correlations(tpc, z_locs, x_locs, y_locs, out_dir)

    # ── 4. Temporal spectra ───────────────────────────────────────────────
    print("\n[4/5] Temporal spectra …")
    tspec = compute_temporal_spectra(G, stations, args.dt_ctu)
    if not args.skip_plots:
        plot_temporal_spectra(tspec, x_locs, y_locs, out_dir)

    save_results(out_dir, stats, spec_z, tpc, tspec, x_locs, y_locs, z_locs)
    print(f"\n✓  Done → {out_dir.resolve()}")


if __name__ == "__main__":
    run(build_parser().parse_args())
