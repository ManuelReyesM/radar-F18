"""
Long-range pulse-Doppler radar vs. a fighter-class aircraft (F/A-18-like target).

Scenario: the aircraft flies a straight, level pass. It is first seen 500 km
from the radar, passes the radar at a closest ground distance of 100 km and
flies away again. Everything is generic textbook radar theory with assumed,
open-literature-style parameters (see the PARAMETERS block).

What the script produces (folder ./figures):
  00_scheme_pass.png        schematic of the flight path (plan + side view)
  01_snr_pd.png             radar equation: SNR and probability of detection
  02_doppler_aliasing.png   true vs. folded Doppler for each PRF
  03_cpi_fft_cfar.png       one pulse-level CPI: FFT + CA-CFAR detection
  04_ambiguity_resolution.png  range / speed recovered from staggered PRFs

Author: Manuel Reyes Martinez Lax
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------- PARAMETERS
C = 3e8
K_B = 1.380649e-23

# Radar (UHF, large phased array, assumed values)
F0 = 450e6                      # carrier [Hz]
LAM = C / F0                    # wavelength [m]
PT = 90e3                       # peak power [W]
G = 10 ** (28 / 10)             # antenna gain (28 dBi)
TAU = 100e-6                    # pulse length [s]
BW = 1e6                        # chirp bandwidth [Hz] -> 150 m range resolution
NF = 10 ** (3 / 10)             # noise figure 3 dB
LOSS = 10 ** (4 / 10)           # system losses 4 dB
T0 = 290.0                      # reference temperature [K]
PRFS = np.array([250.0, 280.0, 310.0])   # staggered PRFs [Hz]
N_PULSES = 32                   # pulses per CPI (coherent integration)
PFA = 1e-6                      # false alarm probability per cell
H_RADAR = 30.0                  # antenna height [m]

# Target
H_T = 14e3                      # flight altitude [m] (near fighter ceiling)
D_CPA = 100e3                   # closest ground distance to radar [m]
R_START = 500e3                 # first position, slant range [m]
SPEEDS = {"Mach 0.85 (250 m/s)": 250.0, "Mach 1.5 (440 m/s)": 440.0}
RCS = 3.0                       # target RCS [m^2] (typical fighter, assumed)

DT = 5.0                        # time between measurement epochs [s]
OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(7)


# --------------------------------------------------------------------- MODEL
def geometry(v):
    """Straight, level flight. Returns time, along-track x, slant range, radial
    speed (positive = approaching) and ground range."""
    x0 = np.sqrt(R_START**2 - D_CPA**2 - H_T**2)
    t = np.arange(0, 2 * x0 / v, DT)
    x = -x0 + v * t
    r = np.sqrt(x**2 + D_CPA**2 + H_T**2)
    vr = -(x * v) / r
    rg = np.sqrt(x**2 + D_CPA**2)
    return t, x, r, vr, rg


def snr_single_pulse(r, sigma):
    num = PT * G**2 * LAM**2 * sigma * TAU
    den = (4 * np.pi) ** 3 * r**4 * K_B * T0 * NF * LOSS
    return num / den


def snr_cpi(r, sigma):
    return snr_single_pulse(r, sigma) * N_PULSES


def pd_swerling1(snr):
    return PFA ** (1.0 / (1.0 + snr))


def horizon():
    """Radar horizon with the 4/3 earth model [m]."""
    return 4.12e3 * (np.sqrt(H_T) + np.sqrt(H_RADAR))


def wrap(x, period):
    return (x + period / 2) % period - period / 2


def resolve_range(r_meas):
    ru = C / (2 * PRFS)
    cand = r_meas[0] + ru[0] * np.arange(0, 9)
    cand = cand[cand <= 800e3]
    err = np.zeros_like(cand)
    for i in range(1, len(PRFS)):
        err += wrap(cand - r_meas[i], ru[i]) ** 2
    return cand[np.argmin(err)]


def resolve_doppler(fd_meas):
    cand = fd_meas[0] + PRFS[0] * np.arange(-12, 13)
    err = np.zeros_like(cand)
    for i in range(1, len(PRFS)):
        err += wrap(cand - fd_meas[i], PRFS[i]) ** 2
    return cand[np.argmin(err)]


def ca_cfar(power, n_train=8, n_guard=2, pfa=1e-3):
    """Cell-averaging CFAR over a circular spectrum."""
    n = len(power)
    alpha = 2 * n_train * (pfa ** (-1.0 / (2 * n_train)) - 1)
    thr = np.zeros(n)
    for i in range(n):
        idx = [(i + k) % n for k in range(-n_guard - n_train, -n_guard)]
        idx += [(i + k) % n for k in range(n_guard + 1, n_guard + n_train + 1)]
        thr[i] = alpha * power[idx].mean()
    return thr


# ------------------------------------------------------------------ FIGURES
def fig_scheme():
    x0 = np.sqrt(R_START**2 - D_CPA**2 - H_T**2) / 1e3
    d = D_CPA / 1e3
    hor = horizon() / 1e3
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 9.5), gridspec_kw={"height_ratios": [1.25, 1]})

    # plan view
    th = np.linspace(-0.6, 0.6, 200)
    for sgn in (1, -1):
        a1.plot(sgn * R_START / 1e3 * np.cos(th), R_START / 1e3 * np.sin(th), "--", color="tab:blue", lw=1,
                label="500 km design range" if sgn == 1 else None)
        a1.plot(sgn * hor * np.cos(th), hor * np.sin(th), ":", color="tab:red", lw=1.3,
                label=f"radar horizon ({hor:.0f} km, 4/3 earth)" if sgn == 1 else None)
    a1.annotate("", xy=(x0, d), xytext=(-x0, d), arrowprops=dict(arrowstyle="-|>", lw=2.4, color="k"))
    a1.plot([-x0, 0, x0], [d, d, d], "ko", ms=6)
    a1.plot(0, 0, "^", ms=14, color="tab:green")
    a1.text(14, -8, "Radar", color="tab:green", fontsize=11, weight="bold", va="center")
    a1.plot([0, 0], [0, d], "-", color="gray", lw=1)
    a1.text(8, d / 2, f"closest distance {d:.0f} km", fontsize=9, va="center")
    a1.plot([0, -x0], [0, d], "-", color="gray", lw=0.8)
    a1.text(-x0 / 2, d / 2 + 14, "R = 500 km", fontsize=9, ha="center", rotation=-11.5)
    a1.text(-x0, d + 20, "START: first contact", ha="center", fontsize=9)
    a1.text(0, d + 20, "CLOSEST POINT: Doppler = 0", ha="center", fontsize=9)
    a1.text(x0, d + 20, "END: leaves range", ha="center", fontsize=9)
    a1.text(-x0 / 2, d + 45, "approaching: fd > 0", ha="center", color="tab:red", fontsize=10)
    a1.text(x0 / 2, d + 45, "receding: fd < 0", ha="center", color="tab:blue", fontsize=10)
    a1.set_aspect("equal")
    a1.set_xlim(-560, 560)
    a1.set_ylim(-40, 190)
    a1.set_xlabel("along-track x [km]")
    a1.set_ylabel("cross-track y [km]")
    a1.set_title("Plan view: straight, level pass of the aircraft")
    a1.legend(loc="lower center", fontsize=8, ncol=2)
    a1.grid(alpha=0.3)

    # side view: minimum altitude the radar can see at each ground distance (4/3 earth)
    dg = np.linspace(0, 560, 400)
    hmin = np.clip((dg - 4.12 * np.sqrt(H_RADAR)) / 4.12, 0, None) ** 2 / 1e3
    a2.fill_between(dg, 0, hmin, color="#cdb78a", alpha=0.8, label="below radar horizon (not visible)")
    a2.axhline(H_T / 1e3, color="k", lw=2.4, label=f"aircraft altitude {H_T/1e3:.0f} km")
    a2.plot(0, 0, "^", color="tab:green", ms=11)
    a2.axvline(hor, color="tab:red", ls=":", label=f"radar horizon at this altitude: {hor:.0f} km")
    a2.axvline(R_START / 1e3, color="tab:blue", ls="--", lw=1, label="500 km design range")
    a2.set_xlim(0, 560); a2.set_ylim(0, 22)
    a2.set_xlabel("ground distance from radar [km]")
    a2.set_ylabel("height above ground [km]")
    a2.set_title("Side view: the aircraft is visible only above the shaded earth-curvature limit")
    a2.legend(fontsize=8, loc="upper left")
    a2.grid(alpha=0.3)
    fig.suptitle("Scenario: UHF long-range radar vs. fighter pass (500 km -> 100 km -> 500 km)")
    fig.tight_layout()
    fig.savefig(OUT / "00_scheme_pass.png", dpi=140)
    plt.close(fig)


def fig_snr_pd():
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    r = np.linspace(50e3, 600e3, 500)
    for s, c in [(1.0, "tab:red"), (3.0, "tab:blue"), (5.0, "tab:green")]:
        ax[0].plot(r / 1e3, 10 * np.log10(snr_cpi(r, s)), color=c, label=f"RCS {s:g} m$^2$")
        ax[1].plot(r / 1e3, pd_swerling1(snr_cpi(r, s)), color=c, label=f"RCS {s:g} m$^2$")
    ax[0].axvline(horizon() / 1e3, color="gray", ls=":", label="radar horizon")
    ax[1].axvline(horizon() / 1e3, color="gray", ls=":")
    ax[0].axhline(13, color="k", lw=0.7, ls="--")
    ax[0].set_xlabel("slant range [km]"); ax[0].set_ylabel("SNR after 32-pulse integration [dB]")
    ax[0].set_title("Radar equation"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("slant range [km]"); ax[1].set_ylabel("Pd (Swerling I, Pfa = 1e-6)")
    ax[1].set_title("Probability of detection"); ax[1].grid(alpha=0.3); ax[1].legend()
    fig.tight_layout(); fig.savefig(OUT / "01_snr_pd.png", dpi=140); plt.close(fig)


def fig_doppler():
    fig, axs = plt.subplots(len(SPEEDS), 1, figsize=(11, 7), sharex=True)
    for ax, (name, v) in zip(axs, SPEEDS.items()):
        t, _, r, vr, _ = geometry(v)
        fd = 2 * vr / LAM
        ax.plot(t / 60, fd, "k", lw=2, label="true Doppler")
        for p, c in zip(PRFS, ["tab:blue", "tab:orange", "tab:green"]):
            ax.plot(t / 60, wrap(fd, p), ".", ms=2.5, color=c, label=f"folded, PRF {p:g} Hz")
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_ylabel("Doppler [Hz]"); ax.set_title(name); ax.grid(alpha=0.3)
        ax.legend(ncol=4, fontsize=8, loc="upper right")
    axs[-1].set_xlabel("time [min]")
    fig.suptitle("Doppler ambiguity: every PRF sees a different aliased copy of the true Doppler")
    fig.tight_layout(); fig.savefig(OUT / "02_doppler_aliasing.png", dpi=140); plt.close(fig)


def fig_cpi(v=250.0):
    """Pulse-level demo of one CPI, as in the lab sessions: slow-time FFT."""
    t, _, r, vr, _ = geometry(v)
    prf = PRFS[1]
    # inbound epoch with R in 150-450 km whose folded Doppler falls mid-spectrum
    inbound = np.where((r < 450e3) & (r > 150e3) & (vr > 0))[0]
    k = inbound[np.argmin(abs(wrap(2 * vr[inbound] / LAM, prf) - 30.0))]
    fd = 2 * vr[k] / LAM
    n = np.arange(N_PULSES)
    snr = snr_cpi(r[k], RCS)
    amp = np.sqrt(snr / N_PULSES)
    noise = (rng.standard_normal(N_PULSES) + 1j * rng.standard_normal(N_PULSES)) / np.sqrt(2)
    s = amp * np.exp(2j * np.pi * fd * n / prf) + noise
    nfft = 256
    spec = np.fft.fftshift(np.fft.fft(s * np.hamming(N_PULSES), nfft))
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / prf))
    p = abs(spec) ** 2
    thr = ca_cfar(p, n_train=40, n_guard=20)
    fig, ax = plt.subplots(2, 1, figsize=(11, 7))
    ax[0].plot(n / prf * 1e3, s.real, ".-", label="I"); ax[0].plot(n / prf * 1e3, s.imag, ".-", label="Q")
    ax[0].set_xlabel("slow time [ms]"); ax[0].set_ylabel("amplitude")
    ax[0].set_title(f"Range gate of the target, R = {r[k]/1e3:.0f} km, SNR = {10*np.log10(snr):.1f} dB")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(f, 10 * np.log10(p / p.max()), color="tab:blue", label="spectrum")
    ax[1].plot(f, 10 * np.log10(thr / p.max()), color="tab:red", label="CA-CFAR threshold")
    det = p > thr
    ax[1].plot(f[det], 10 * np.log10(p[det] / p.max()), "o", color="tab:green", ms=3, label="CFAR detections")
    ax[1].axvline(wrap(fd, prf), color="k", ls="--", label=f"true Doppler {fd:.0f} Hz, folded to {wrap(fd, prf):.0f} Hz")
    ax[1].set_xlabel("Doppler [Hz]"); ax[1].set_ylabel("dB"); ax[1].set_ylim(-45, 3)
    ax[1].set_title(f"FFT over {N_PULSES} pulses, PRF {prf:g} Hz (unambiguous Doppler only +/-{prf/2:g} Hz)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(OUT / "03_cpi_fft_cfar.png", dpi=140); plt.close(fig)


def track(v):
    """Measure at every epoch with all PRFs, then resolve ambiguities."""
    t, _, r, vr, rg = geometry(v)
    fd_true = 2 * vr / LAM
    snr = snr_cpi(r, RCS)
    pd = pd_swerling1(snr)
    visible = rg <= horizon()
    det = (rng.random(len(t)) < pd) & visible
    ru = C / (2 * PRFS)
    r_res = np.full(len(t), np.nan)
    v_res = np.full(len(t), np.nan)
    r_single = np.full(len(t), np.nan)
    v_single = np.full(len(t), np.nan)
    for i in np.where(det)[0]:
        sr = (C / (2 * BW)) / np.sqrt(2 * snr[i])
        rm = np.array([(r[i] + sr * rng.standard_normal()) % ru[j] for j in range(len(PRFS))])
        fm = np.array([wrap(fd_true[i] + (PRFS[j] / N_PULSES) / np.sqrt(2 * snr[i]) * rng.standard_normal(), PRFS[j])
                       for j in range(len(PRFS))])
        r_res[i] = resolve_range(rm)
        v_res[i] = resolve_doppler(fm) * LAM / 2
        r_single[i] = rm[1]
        v_single[i] = fm[1] * LAM / 2
    return dict(t=t, r=r, vr=vr, det=det, r_res=r_res, v_res=v_res, r_single=r_single, v_single=v_single, pd=pd)


def fig_resolution():
    fig, axs = plt.subplots(2, len(SPEEDS), figsize=(14, 7), sharex="col")
    for col, (name, v) in enumerate(SPEEDS.items()):
        d = track(v)
        tm = d["t"] / 60
        axs[0, col].plot(tm, d["r"] / 1e3, "k", lw=2, label="true")
        axs[0, col].plot(tm, d["r_single"] / 1e3, ".", ms=3, color="tab:gray", label=f"single PRF {PRFS[1]:g} Hz")
        axs[0, col].plot(tm, d["r_res"] / 1e3, ".", ms=4, color="tab:red", label="staggered PRFs, resolved")
        axs[0, col].set_ylabel("range [km]"); axs[0, col].set_title(name); axs[0, col].legend(fontsize=8)
        axs[1, col].plot(tm, d["vr"], "k", lw=2, label="true")
        axs[1, col].plot(tm, d["v_single"], ".", ms=3, color="tab:gray", label="single PRF")
        axs[1, col].plot(tm, d["v_res"], ".", ms=4, color="tab:red", label="staggered PRFs, resolved")
        axs[1, col].set_ylabel("radial speed [m/s]"); axs[1, col].set_xlabel("time [min]")
        axs[1, col].legend(fontsize=8)
        for a in axs[:, col]:
            a.grid(alpha=0.3)
        ok = ~np.isnan(d["r_res"])
        e_r = np.abs(d["r_res"][ok] - d["r"][ok])
        e_v = np.abs(d["v_res"][ok] - d["vr"][ok])
        print(f"{name}: {ok.sum()} detections of {len(d['t'])} epochs | "
              f"resolved range error rms {np.sqrt(np.mean(e_r**2)):.0f} m, speed error rms {np.sqrt(np.mean(e_v**2)):.2f} m/s | "
              f"gross errors range {np.mean(e_r > 1e3)*100:.1f} %, speed {np.mean(e_v > 5)*100:.1f} %")
    fig.suptitle("Ambiguity resolution with three staggered PRFs")
    fig.tight_layout(); fig.savefig(OUT / "04_ambiguity_resolution.png", dpi=140); plt.close(fig)


def print_tables():
    print("\n--- Design numbers ---")
    print(f"lambda = {LAM:.3f} m, unambiguous range per PRF [km]:",
          np.round(C / (2 * PRFS) / 1e3, 0))
    print(f"single-PRF unambiguous speed [m/s] (lambda*PRF/4):", np.round(LAM * PRFS / 4, 1))
    print(f"radar horizon vs. aircraft at {H_T/1e3:.0f} km: {horizon()/1e3:.0f} km")
    for v in SPEEDS.values():
        print(f"v = {v:.0f} m/s -> wavelength for ambiguity-free R=500 km with a single PRF: "
              f"{8 * R_START * v / C:.2f} m (VHF/HF region)")
    for s in (1.0, 3.0, 5.0):
        r = np.linspace(50e3, 1000e3, 40000)
        ok = r[pd_swerling1(snr_cpi(r, s)) >= 0.5]
        print(f"RCS {s:g} m^2: Pd >= 0.5 up to {ok.max()/1e3:.0f} km")


if __name__ == "__main__":
    fig_scheme()
    fig_snr_pd()
    fig_doppler()
    fig_cpi()
    fig_resolution()
    print_tables()
    print(f"\nFigures written to {OUT}")
