#!/usr/bin/env python3
"""
Inundation time-series figures for the Lake Chilwa manuscript.

Reads the pre-computed annual and monthly CSVs in 03.outputs and writes
publication figures F3-F6 to 04.images. No Earth Engine dependency: this is a
post-processing step over exported summaries, so it runs locally without rgee.

Figures
  F3  Sensor availability by year, 1984-2024 (optical Landsat + Sentinel-1 SAR)
  F4  Multi-decadal inundation hydrograph (open water + flooded vegetation)
  F5  Spectral-index behaviour and agreement with the SMA water record
  F6  SAR dynamics (Sentinel-1 seasonality and interannual; PALSAR L-band)

Paths are resolved relative to this script's location.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "03.outputs"))
IMG = os.path.normpath(os.path.join(HERE, "..", "04.images"))
os.makedirs(IMG, exist_ok=True)

# ----------------------------------------------------------------------------
# Shared publication theme
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.8,
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "legend.frameon": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

C = {
    "water": "#2166ac", "fveg": "#5aae61",
    "L5 TM": "#4393c3", "L7 ETM+": "#2166ac", "L8 OLI": "#053061", "S1 C-band": "#b2182b",
    "NDWI": "#4393c3", "MNDWI": "#2166ac", "AWEIsh": "#053061", "WRI": "#d6604d", "NDPI": "#f4a582",
    "VV": "#762a83", "VH": "#9970ab",
    "HH": "#b35806", "HV": "#e08214", "HHHV": "#542788",
    "grid": "#e6e6e6", "accent": "#b2182b",
}
YRS = (1983.5, 2024.5)


def _grid(ax, axis="y"):
    ax.grid(axis=axis, color=C["grid"], linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def _panel_tag(ax, tag):
    ax.text(-0.02, 1.03, tag, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="right")


def savefig(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(IMG, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
sma = pd.read_csv(os.path.join(OUT, "sma_annual_fractions.csv"))
idx = pd.read_csv(os.path.join(OUT, "landsat_annual_indices.csv"))
s1 = pd.read_csv(os.path.join(OUT, "s1_monthly_timeseries.csv"))
pal = pd.read_csv(os.path.join(OUT, "palsar_lband_annual.csv"))
sen = pd.read_csv(os.path.join(OUT, "sensor_availability_by_year.csv"))
sen = sen[sen["year"] <= 2024].copy()
# The audit reports each sensor twice, once unscreened and once cloud-screened,
# and Sentinel-1 only unscreened. Take the unscreened rows so the panel measures
# acquisition availability on one basis and no scene is counted twice.
if "screen" in sen.columns:
    sen = sen[sen["screen"] == "all"].copy()

print("\n================ SUMMARY NUMBERS ================")

# ----------------------------------------------------------------------------
# F4  Inundation hydrograph
# ----------------------------------------------------------------------------
def fig4():
    fig, (ax, axn) = plt.subplots(
        2, 1, figsize=(9, 5.2), sharex=True,
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.08})

    ax.plot(sma.year, sma.water_frac, "-o", color=C["water"], ms=3.4, lw=1.8,
            label="Open water", zorder=3)
    ax.plot(sma.year, sma.flooded_veg_frac, "-s", color=C["fveg"], ms=3.0, lw=1.4,
            label="Emergent / flooded vegetation", zorder=3)

    wmean = sma.water_frac.mean()
    ax.axhline(wmean, color=C["water"], lw=0.8, ls=":", alpha=0.7)
    ax.text(2024.3, wmean, f" mean {wmean:.3f}", color=C["water"], va="center",
            ha="left", fontsize=8)

    # mark principal recession troughs: local minima below a tenth of the basin
    w = sma.reset_index(drop=True)
    v = w.water_frac.values
    li = [i for i in range(1, len(v) - 1)
          if v[i] < v[i - 1] and v[i] < v[i + 1] and v[i] < 0.09]
    troughs = w.iloc[li]
    ax.scatter(troughs.year, troughs.water_frac, s=46, facecolors="none",
               edgecolors=C["accent"], linewidths=1.2, zorder=4)
    for _, r in troughs.iterrows():
        ax.annotate(f"{int(r.year)}", (r.year, r.water_frac),
                    textcoords="offset points", xytext=(0, -13),
                    ha="center", fontsize=8, color=C["accent"])
    hi = w.loc[w.water_frac.idxmax()]
    ax.annotate("record high\n2023-24", (hi.year, hi.water_frac),
                textcoords="offset points", xytext=(-8, 2), ha="right",
                va="bottom", fontsize=8, color=C["water"])

    _grid(ax)
    ax.set_ylabel("Basin fraction")
    ax.set_ylim(0, max(sma.water_frac.max(), sma.flooded_veg_frac.max()) * 1.18)
    ax.legend(loc="upper left", ncol=2)
    ax.set_title("Reconstructed inundation of the Lake Chilwa basin, 1984-2024",
                 loc="left")
    _panel_tag(ax, "a")

    axn.bar(sma.year, sma.n_scenes, color="#9e9e9e", width=0.7, zorder=2)
    _grid(axn)
    axn.set_ylabel("Scenes", fontsize=8.5)
    axn.set_xlabel("Year")
    axn.set_xlim(*YRS)
    axn.xaxis.set_major_locator(MultipleLocator(5))
    _panel_tag(axn, "b")
    savefig(fig, "fig04_inundation_hydrograph")

    wmin = sma.loc[sma.water_frac.idxmin()]
    wmax = sma.loc[sma.water_frac.idxmax()]
    print(f"F4 water_frac: mean={wmean:.3f}, min={wmin.water_frac:.3f} ({int(wmin.year)}), "
          f"max={wmax.water_frac:.3f} ({int(wmax.year)})")
    print(f"   principal troughs (local min < 0.09): {[int(y) for y in troughs.year]}")
    print(f"   flooded_veg: mean={sma.flooded_veg_frac.mean():.3f}, "
          f"min={sma.flooded_veg_frac.min():.3f}, max={sma.flooded_veg_frac.max():.3f}")


# ----------------------------------------------------------------------------
# F5  Spectral-index behaviour and agreement with SMA
# ----------------------------------------------------------------------------
def fig5():
    m = idx.merge(sma[["year", "water_frac"]], on="year")
    cols = ["NDWI", "MNDWI", "AWEIsh", "WRI", "NDPI"]

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(10.5, 4.4), gridspec_kw={"width_ratios": [2.3, 1]})

    for c in cols:
        z = (m[c] - m[c].mean()) / m[c].std()
        ax.plot(m.year, z, "-", color=C[c], lw=1.5, label=c)
    zw = (m.water_frac - m.water_frac.mean()) / m.water_frac.std()
    ax.plot(m.year, zw, "--", color="black", lw=1.6, label="SMA open water")
    _grid(ax)
    ax.set_ylabel("Standardised anomaly (z)")
    ax.set_xlabel("Year")
    ax.set_xlim(*YRS)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.legend(loc="upper left", ncol=3)
    ax.set_title("Index anomalies against the SMA water record", loc="left")
    _panel_tag(ax, "a")

    cors = {c: np.corrcoef(m[c], m.water_frac)[0, 1] for c in cols}
    order = sorted(cols, key=lambda c: cors[c])
    bars = axb.barh([o for o in order], [cors[o] for o in order],
                    color=[C[o] for o in order], zorder=2)
    _grid(axb, axis="x")
    axb.axvline(0, color="#444444", lw=0.8)
    axb.set_xlabel("Pearson r with SMA water")
    axb.set_xlim(-1.3, 1.2)
    axb.tick_params(axis="y", pad=6)
    for o in order:
        pos = cors[o] >= 0
        axb.text(cors[o] + (0.03 if pos else -0.03), o, f"{cors[o]:+.2f}",
                 va="center", ha="left" if pos else "right", fontsize=8.5)
    axb.set_title("Agreement with inundation", loc="left")
    _panel_tag(axb, "b")
    savefig(fig, "fig05_index_performance")

    r_np = np.corrcoef(m.NDPI, m.MNDWI)[0, 1]
    print("F5 index-vs-water Pearson r: " + ", ".join(f"{c}={cors[c]:+.2f}" for c in cols))
    print(f"   NDPI vs MNDWI r={r_np:+.3f}  (redundant: NDPI = -MNDWI by construction)")


# ----------------------------------------------------------------------------
# F6  SAR dynamics
# ----------------------------------------------------------------------------
def fig6():
    fig = plt.figure(figsize=(10.5, 7.4))
    gs = fig.add_gridspec(3, 1, hspace=0.42)

    # (a) Sentinel-1 monthly climatology
    axa = fig.add_subplot(gs[0])
    clim = s1.groupby("month").agg(
        VV=("VV_mean", "mean"), VH=("VH_mean", "mean"),
        wf=("water_frac", "mean"), wf_lo=("water_frac", "min"),
        wf_hi=("water_frac", "max")).reset_index()
    mo = clim.month
    axa.plot(mo, clim.VV, "-o", color=C["VV"], ms=3.5, lw=1.6, label="VV")
    axa.plot(mo, clim.VH, "-o", color=C["VH"], ms=3.5, lw=1.6, label="VH")
    axa.set_ylabel("Backscatter (dB)")
    axa.set_xticks(range(1, 13))
    axa.set_xticklabels(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
    _grid(axa)
    axr = axa.twinx()
    axr.spines["top"].set_visible(False)
    axr.fill_between(mo, clim.wf_lo, clim.wf_hi, color=C["water"], alpha=0.12)
    axr.plot(mo, clim.wf, "-s", color=C["water"], ms=3.2, lw=1.6, label="Water fraction")
    axr.set_ylabel("Water fraction", color=C["water"])
    axr.tick_params(axis="y", colors=C["water"])
    h1, l1 = axa.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    axa.legend(h1 + h2, l1 + l2, loc="lower left", ncol=3)
    axa.set_title("Sentinel-1 seasonal climatology (2015-2024)", loc="left")
    _panel_tag(axa, "a")

    # (b) Sentinel-1 interannual monthly water fraction
    axb = fig.add_subplot(gs[1])
    s1s = s1.assign(t=pd.to_datetime(s1.date)).sort_values("t")
    axb.plot(s1s.t, s1s.water_frac, "-", color=C["water"], lw=1.3)
    axb.scatter(s1s.t, s1s.water_frac, s=9, color=C["water"], zorder=3)
    _grid(axb)
    axb.set_ylabel("Water fraction")
    axb.set_title("Sentinel-1 monthly water fraction", loc="left")
    _panel_tag(axb, "b")

    # (c) PALSAR L-band annual (average duplicate acquisitions per year)
    axc = fig.add_subplot(gs[2])
    pa = pal.groupby("year").mean().reindex(range(pal.year.min(), pal.year.max() + 1))
    axc.plot(pa.index, pa.L_HH, "-o", color=C["HH"], ms=3.5, lw=1.5, label="HH")
    axc.plot(pa.index, pa.L_HV, "-o", color=C["HV"], ms=3.5, lw=1.5, label="HV")
    _grid(axc)
    axc.set_ylabel("Backscatter (dB)")
    axc.set_xlabel("Year")
    axc.set_xticks(range(2008, 2025, 2))
    axc.set_ylim(top=axc.get_ylim()[1] + 1.4)
    axcr = axc.twinx()
    axcr.spines["top"].set_visible(False)
    axcr.plot(pa.index, pa.L_HH_HV, "-^", color=C["HHHV"], ms=4, lw=1.6,
              label="HH - HV (double bounce)")
    axcr.set_ylabel("HH - HV (dB)", color=C["HHHV"])
    axcr.tick_params(axis="y", colors=C["HHHV"])
    hc1, lc1 = axc.get_legend_handles_labels()
    hc2, lc2 = axcr.get_legend_handles_labels()
    axc.legend(hc1 + hc2, lc1 + lc2, loc="upper center", ncol=3)
    axc.set_title("PALSAR L-band annual backscatter (flooded-vegetation sensitivity)",
                  loc="left")
    _panel_tag(axc, "c")
    savefig(fig, "fig06_sar_seasonality")

    hi = clim.loc[clim.wf.idxmax()]
    lo = clim.loc[clim.wf.idxmin()]
    print(f"F6 S1 climatology water fraction: high month={int(hi.month)} ({hi.wf:.3f}), "
          f"low month={int(lo.month)} ({lo.wf:.3f})")
    print(f"   PALSAR HH-HV mean={pal.L_HH_HV.mean():.2f} dB "
          f"({pal.year.min()}-{pal.year.max()}, gap 2011-2014)")


# ----------------------------------------------------------------------------
# F3  Sensor availability
# ----------------------------------------------------------------------------
def fig3():
    order = ["L5 TM", "L7 ETM+", "L8 OLI", "S1 C-band"]
    piv = (sen.groupby(["year", "sensor"])["n"].sum().unstack("sensor")
           .reindex(columns=order).reindex(range(1984, 2025)).fillna(0))
    fig, ax = plt.subplots(figsize=(9.6, 4.2))
    bottom = np.zeros(len(piv))
    for s in order:
        ax.bar(piv.index, piv[s], bottom=bottom, color=C[s], width=0.8,
               label=s, zorder=2)
        bottom += piv[s].values
    _grid(ax)
    ax.set_ylabel("Scenes per year")
    ax.set_xlabel("Year")
    ax.set_xlim(*YRS)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.legend(loc="upper left", ncol=4)
    ax.set_title("Optical and radar scene availability over the Lake Chilwa basin",
                 loc="left")
    savefig(fig, "fig03_sensor_availability")

    tot = piv.sum(axis=1)
    print(f"F3 scenes/yr: 1984-1999 mean={tot.loc[1984:1999].mean():.0f}, "
          f"2015-2024 mean={tot.loc[2015:2024].mean():.0f}")
    print("   sensors:", ", ".join(order))


if __name__ == "__main__":
    print("Generating figures ->", IMG)
    fig4(); fig5(); fig6(); fig3()
    print("================================================\n")
