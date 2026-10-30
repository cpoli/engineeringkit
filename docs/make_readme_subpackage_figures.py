"""Regenerate the README's per-subpackage teaser figures.

    python docs/make_readme_subpackage_figures.py                  # all subpackages
    python docs/make_readme_subpackage_figures.py materials aerodynamics  # just these

Writes docs/source/_static/images/readme_<subpackage>.png, three panels each,
at the same size as readme_hero.png. README.md embeds them by their
raw.githubusercontent.com URLs so they also render on PyPI.
``hero`` writes readme_hero.png, the banner on the docs landing page.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).parent / "source" / "_static" / "images"


def structural_engineering(axes):
    from engineeringkit.structural_engineering import Beam, PlanarFrame, method_of_joints
    from engineeringkit.structural_engineering.visualizers import plot_frame, plot_truss

    ax1, ax2, ax3 = axes
    panels, a, depth, P = 6, 4.0, 4.0, 50e3
    nodes = np.array([(a * k, 0.0) for k in range(panels + 1)] + [(a * k, depth) for k in range(1, panels)])
    top = {k: panels + k for k in range(1, panels)}
    members = [(k, k + 1) for k in range(panels)] + [(top[k], top[k + 1]) for k in range(1, panels - 1)]
    members += [(0, top[1]), (panels, top[panels - 1])] + [(k, top[k]) for k in range(1, panels)]
    members += [(top[k], k + 1) for k in range(1, panels // 2)] + [(top[k], k - 1) for k in range(panels // 2 + 1, panels)]
    loads = np.zeros(2 * len(nodes))
    loads[2 * np.arange(1, panels) + 1] = -P
    result = method_of_joints(nodes, members, loads, [0, 1, 2 * panels + 1])
    plot_truss(nodes, members, forces=result.member_forces, ax=ax1)
    ax1.set_title("Pratt truss: tension red, compression blue")

    beam = Beam(12.0, 200e9, 8e-5, ((0.0, "pin"), (5.0, "roller"), (12.0, "roller")))
    r = beam.solve(point_loads=((9.0, 30e3),), distributed_loads=((0.0, 12.0, 10e3),))
    ax2.fill_between(r.x, r.moment / 1e3, alpha=0.3, color="C0")
    ax2.plot(r.x, r.moment / 1e3, color="C0")
    ax2.axhline(0, color="k", lw=0.6)
    ax2.plot([0, 5, 12], [0, 0, 0], "k^", ms=9)
    ax2.set(xlabel="x [m]", ylabel="bending moment [kN m]", title="Two-span continuous beam (Macaulay)")
    twin = ax2.twinx()
    twin.plot(r.x, 1e3 * r.deflection, color="C3")
    twin.set_ylabel("deflection [mm]", color="C3")
    ax2.set_ylim(np.array([-1.1, 1.1]) * np.abs(r.moment / 1e3).max())
    twin.set_ylim(np.array([-1.1, 1.1]) * np.abs(1e3 * r.deflection).max())

    h, L = 4.0, 6.0
    frame = PlanarFrame([[0, 0], [0, h], [L, h], [L, 0]], ((0, 1), (1, 2), (2, 3)), 200e9, 1e-2, 2e-4)
    frame_loads = np.zeros(12)
    frame_loads[3] = 40e3
    plot_frame(frame.solve(frame_loads, [0, 1, 2, 9, 10, 11]), scale=60, ax=ax3)
    ax3.set_title("Portal frame sway (stiffness method)")


def mechanical_engineering(axes):
    from engineeringkit.mechanical_engineering import fourbar_kinematics, normal_modes, transmissibility
    from engineeringkit.mechanical_engineering.visualizers import plot_frequency_response, plot_linkage, plot_mode_shapes

    ax1, ax2, ax3 = axes
    r = np.logspace(-1, 1, 600)
    for zeta in (0.05, 0.1, 0.25, 0.5, 1.0):
        plot_frequency_response(r, transmissibility(r, zeta), ax=ax1, label=rf"$\zeta$ = {zeta}")
    ax1.axvline(np.sqrt(2), color="gray", ls=":")
    ax1.set_yscale("log")
    ax1.set_title(r"Vibration isolation: isolates only above $\sqrt{2}\,\omega_n$")
    ax1.legend(fontsize=8)

    n = 5
    M = np.eye(n)
    K = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    K[-1, -1] = 1.0
    modes = normal_modes(M, K)
    plot_mode_shapes(modes.mode_shapes[:3], modes.frequencies[:3], ax=ax2)
    ax2.set_title("Normal modes of a five-mass chain")

    a, b, c, d = 1.0, 3.0, 2.5, 3.0
    theta2 = np.linspace(0, 2 * np.pi, 361)
    motion = fourbar_kinematics(a, b, c, d, theta2, 1.0)
    A = a * np.exp(1j * theta2)
    B = d + c * np.exp(1j * motion.theta4)
    coupler = A + 0.6 * (B - A) + 0.8j * (B - A) / abs(B - A)
    ax3.plot(coupler.real, coupler.imag, color="C3", label="coupler curve")
    for k in (0, 120, 240):
        plot_linkage([[0, 0], [A[k].real, A[k].imag], [B[k].real, B[k].imag], [d, 0]], ax=ax3, color="C0", alpha=0.4)
    ax3.set_aspect("equal")
    ax3.set_title("Crank-rocker four-bar and its coupler curve")
    ax3.legend(fontsize=8, loc="lower right")


def thermodynamics(axes):
    from engineeringkit.thermodynamics import (
        brayton_efficiency,
        carnot_efficiency,
        diesel_efficiency,
        isentropic_pressure_volume,
        otto_efficiency,
        steam_properties,
    )
    from engineeringkit.thermodynamics.visualizers import plot_pv_diagram, plot_ts_diagram

    ax1, ax2, ax3 = axes
    s4 = steam_properties(10e3, x=0.0)
    s1 = steam_properties(10e6, s=s4.s)
    s2 = steam_properties(10e6, 773.15)
    s3 = steam_properties(10e3, s=s2.s)
    T_boil = steam_properties(10e6, x=0.0).T
    path = (
        [steam_properties(10e6, T) for T in np.linspace(s1.T, T_boil - 0.01, 60)]
        + [steam_properties(10e6, x=x) for x in np.linspace(0, 1, 20)]
        + [steam_properties(10e6, T) for T in np.linspace(T_boil + 0.01, s2.T, 60)]
    )
    s_cycle = [s4.s] + [st.s for st in path] + [s3.s, s4.s]
    T_cycle = [s4.T] + [st.T for st in path] + [s3.T, s4.T]
    plot_ts_diagram(s_cycle, T_cycle, ax=ax1, label="ideal Rankine cycle")
    ax1.set_title("Rankine cycle on the steam dome (IAPWS-IF97)")
    ax1.legend(fontsize=8, loc="upper left")

    r, alpha, rc, g, p1 = 16.0, 1.5, 1.4, 1.4, 100e3
    v2, v3 = 1.0 / r, rc / r
    compression, expansion = np.linspace(1.0, v2, 40), np.linspace(v3, 1.0, 40)
    p3 = alpha * p1 * r**g
    v = np.r_[compression, v2, v3, expansion, 1.0]
    p = np.r_[isentropic_pressure_volume(p1, 1.0, compression, g), p3, p3, isentropic_pressure_volume(p3, v3, expansion, g), p1]
    plot_pv_diagram(v, p, ax=ax2, label="dual cycle")
    ax2.set(xlabel="volume / initial volume", title="Dual cycle, compression ratio 16")

    ratio = np.linspace(2, 24, 200)
    ax3.plot(ratio, otto_efficiency(ratio), label="Otto")
    ax3.plot(ratio, diesel_efficiency(ratio, 2.0), label="Diesel, cut-off 2")
    ax3.plot(ratio, brayton_efficiency(ratio), label="Brayton (pressure ratio)")
    ax3.axhline(carnot_efficiency(1800.0, 300.0), color="gray", ls=":", label="Carnot, 1800 K / 300 K")
    ax3.set(xlabel="compression or pressure ratio", ylabel="efficiency", ylim=(0, 1), title="Air-standard cycle efficiencies")
    ax3.legend(fontsize=8, loc="lower right")


def fluid_engineering(axes):
    from engineeringkit.fluid_engineering import critical_depth, pump_curve_head, pump_operating_point, specific_energy, system_curve_head
    from engineeringkit.fluid_engineering.visualizers import plot_moody_diagram

    ax1, ax2, ax3 = axes
    plot_moody_diagram(ax=ax1)
    ax1.set_title("Moody diagram (Colebrook)")
    ax1.legend(fontsize=7, loc="lower left", ncol=2)

    Q = np.linspace(0, 0.12, 200)
    H0, a, Hs, k = 40.0, 3000.0, 10.0, 2000.0
    ax2.plot(1e3 * Q, system_curve_head(Q, Hs, k), "k", lw=2, label="system")
    for series, parallel, label in ((1, 1, "one pump"), (2, 1, "two in series"), (1, 2, "two in parallel")):
        Qp = Q / parallel
        ax2.plot(1e3 * Q, series * pump_curve_head(Qp, H0, a), label=label)
        point = pump_operating_point(H0, a, Hs, k, series=series, parallel=parallel)
        ax2.plot(1e3 * point.flow, point.head, "ko", ms=5)
    ax2.set(xlabel="flow rate [L/s]", ylabel="head [m]", ylim=(0, 90), title="Pump operating points")
    ax2.legend(fontsize=8)

    y = np.linspace(0.15, 3.0, 400)
    for q in (1.0, 2.0, 4.0):
        ax3.plot(specific_energy(y, q), y, label=f"q = {q:g} m²/s")
        yc = critical_depth(q)
        ax3.plot(specific_energy(yc, q), yc, "ko", ms=4)
    ax3.plot([0, 3.0], [0, 3.0], "k:", lw=0.8)
    ax3.set(xlabel="specific energy E [m]", ylabel="depth y [m]", xlim=(0, 3.5), ylim=(0, 3), title="Specific energy curves")
    ax3.legend(fontsize=8)


def heat_transfer(axes):
    from engineeringkit.heat_transfer import conduction_2d, exchanger_effectiveness_counterflow, exchanger_effectiveness_parallel
    from engineeringkit.heat_transfer.visualizers import plot_blackbody_spectra

    ax1, ax2, ax3 = axes
    plate = conduction_2d(2.0, 1.0, 81, 41, 300.0, 300.0, 300.0, 400.0, k=15.0)
    plate.plot(ax=ax1)
    ax1.set_title("Steady 2-D conduction, top edge at 400 K")

    plot_blackbody_spectra([3000.0, 4000.0, 5000.0, 5800.0], ax=ax2)
    ax2.set_title("Planck spectra with Wien's displacement law")

    ntu = np.linspace(0, 5, 200)
    for Cr, color in ((0.0, "k"), (0.5, "C0"), (1.0, "C3")):
        ax3.plot(ntu, exchanger_effectiveness_counterflow(ntu, Cr), color=color, label=f"counterflow, $C_r$ = {Cr:g}")
        if Cr:
            ax3.plot(ntu, exchanger_effectiveness_parallel(ntu, Cr), "--", color=color, label=f"parallel, $C_r$ = {Cr:g}")
    ax3.set(xlabel="NTU", ylabel="effectiveness", ylim=(0, 1.02), title="Heat exchanger effectiveness-NTU")
    ax3.legend(fontsize=8, loc="lower right")


def aerodynamics(axes):
    from engineeringkit.aerodynamics import lifting_line, naca4_airfoil, panel_method
    from engineeringkit.aerodynamics.visualizers import plot_pressure_coefficient, plot_spanwise_loading, plot_theta_beta_mach

    ax1, ax2, ax3 = axes
    plot_pressure_coefficient(panel_method(naca4_airfoil("2412", 121), np.radians(4.0)), ax=ax1)
    ax1.set_title(r"NACA 2412 at 4°: panel-method $-C_p$")

    b = 10.0
    wing = lifting_line(b, lambda y: 1.6 * (1 - 0.7 * np.abs(2 * y / b)), np.radians(5.0))
    plot_spanwise_loading(wing, ax=ax2)
    ax2.set_title(f"Lifting line: tapered wing, e = {wing.span_efficiency:.3f}")

    plot_theta_beta_mach([1.5, 2.0, 3.0, 5.0, 10.0], ax=ax3)
    ax3.set(xlim=(0, 45), ylim=(0, 90), title="Oblique shocks: θ-β-M")


def electrical_engineering(axes):
    from engineeringkit.electrical_engineering import (
        induction_motor_torque,
        input_impedance,
        reflection_coefficient,
        rlc_damping,
        rlc_step_voltage,
    )
    from engineeringkit.electrical_engineering.visualizers import plot_smith_chart

    ax1, ax2, ax3 = axes
    L, C = 1e-3, 1e-6
    R_critical = 2 * np.sqrt(L / C)
    t = np.linspace(0, 1.2e-3, 2001)
    for R in (4 * R_critical, R_critical, 0.1 * R_critical):
        ax1.plot(1e3 * t, rlc_step_voltage(1.0, R, L, C, t), label=rlc_damping(R, L, C)[2])
    ax1.set(xlabel="time [ms]", ylabel="capacitor voltage [V]", title="Series RLC step response")
    ax1.legend(fontsize=8)

    z0, zl = 50.0, 25.0 - 30.0j
    locus = reflection_coefficient(input_impedance(zl, z0, np.linspace(0, np.pi, 361)), z0)
    plot_smith_chart(locus, ax=ax2, label="along the line")
    gamma_l = reflection_coefficient(zl, z0)
    ax2.plot(gamma_l.real, gamma_l.imag, "ks", label="load")
    ax2.legend(loc="upper left", fontsize=8)
    ax2.set_title("Smith chart: a load seen along a line")

    v_phase, ws = 400 / np.sqrt(3), 2 * np.pi * 50.0 / 2
    slip = np.linspace(1e-4, 1, 1000)
    speed = (1 - slip) * ws * 60 / (2 * np.pi)
    for r2 in (0.25, 0.6, 1.5):
        ax3.plot(speed, induction_motor_torque(v_phase, slip, ws, 0.3, 0.6, r2, 0.6, 20.0), label=f"$R_2$ = {r2} Ω")
    ax3.set(xlabel="speed [rpm]", ylabel="torque [N m]", title="Induction motor torque-speed curves")
    ax3.legend(fontsize=8)


def control_systems(axes):
    from engineeringkit.control_systems import PID, transfer_function
    from engineeringkit.control_systems.visualizers import plot_nyquist, plot_response, plot_root_locus

    ax1, ax2, ax3 = axes
    plant = transfer_function([1], [1, 3, 3, 1])
    for kp, ki, kd in ((1.0, 0.0, 0.0), (2.0, 0.5, 0.0), (2.5, 1.0, 1.5)):
        r = PID(kp, ki, kd).simulate(plant, T=np.linspace(0, 20, 1000))
        plot_response(r.t, r.y, ax=ax1, label=f"$K_p$={kp}, $K_i$={ki}, $K_d$={kd}")
    ax1.axhline(1, color="gray", lw=0.6)
    ax1.set_title(r"PID control of $1/(s+1)^3$")
    ax1.legend(fontsize=8, loc="lower right")

    num, den = [1.0], [1.0, 6.0, 8.0, 0.0]
    plot_root_locus(num, den, np.concatenate((np.linspace(0, 5, 400), np.linspace(5, 200, 800))), ax=ax2)
    ax2.plot([0, 0], [-np.sqrt(8), np.sqrt(8)], "C3o", label="$K = 48$: marginal")
    ax2.set(xlim=(-8, 3), ylim=(-6, 6), title=r"Root locus of $K/(s(s+2)(s+4))$")
    ax2.legend(fontsize=8)

    plot_nyquist(transfer_function([2], [1, 3, 2, 0]), np.logspace(-1, 2, 2000), ax=ax3)
    ax3.set(xlim=(-2, 0.5), ylim=(-2, 2), title=r"Nyquist plot of $2/(s(s+1)(s+2))$")


def signal_processing(axes):
    from scipy import signal

    from engineeringkit.signal_processing import butterworth_lowpass, chebyshev1_lowpass, elliptic_lowpass, fourier_partial_sum, stft
    from engineeringkit.signal_processing.visualizers import plot_magnitude_response

    ax1, ax2, ax3 = axes
    t = np.linspace(-0.25, 0.25, 4001)
    for K in (5, 25, 125):
        k = np.arange(1, K + 1)
        ax1.plot(t, fourier_partial_sum(t, 1.0, np.zeros(K), np.where(k % 2 == 1, 4 / (np.pi * k), 0.0)), lw=0.8, label=f"K = {K}")
    ax1.plot(t, np.sign(t), "k", lw=1.5, label="square wave")
    ax1.set(xlabel="t [s]", title="Gibbs phenomenon: a 9% overshoot that stays")
    ax1.legend(fontsize=8)

    fs, fc = 1000.0, 100.0
    designs = {
        "Butterworth": butterworth_lowpass(6, fc, fs),
        "Chebyshev I": chebyshev1_lowpass(6, 1.0, fc, fs),
        "elliptic": elliptic_lowpass(6, 1.0, 60.0, fc, fs),
    }
    for label, design in designs.items():
        plot_magnitude_response(*design, fs=fs, ax=ax2, label=label)
    ax2.set(ylim=(-90, 5), title="Sixth-order low-pass designs")
    ax2.legend(fontsize=8)

    ts = np.arange(int(2 * fs)) / fs
    chirp = signal.chirp(ts, 50, ts[-1], 400)
    chirp[int(1.5 * fs)] += 20.0
    stft(chirp, fs, 128, hop=8, window=("gaussian", 128 / 6)).plot(ax=ax3)
    ax3.set_title("STFT of a chirp with a click")


def communications(axes):
    from engineeringkit.communications import (
        ber_bpsk_awgn,
        ber_bpsk_rayleigh,
        bpsk,
        qam_constellation,
        root_raised_cosine_pulse,
        ser_mqam_awgn,
    )
    from engineeringkit.communications.visualizers import plot_ber_curve, plot_constellation, plot_eye_diagram

    ax1, ax2, ax3 = axes
    rng = np.random.default_rng(1948)
    points = qam_constellation(16)
    points = points / np.sqrt(np.mean(np.abs(points) ** 2))
    tx = points[rng.integers(0, 16, 4000)]
    plot_constellation(tx + (rng.normal(0, 0.08, tx.size) + 1j * rng.normal(0, 0.08, tx.size)), ax=ax1)
    ax1.plot(points.real, points.imag, "k+", ms=10)
    ax1.set_title("16-QAM through additive white Gaussian noise")

    sps, span = 16, 8
    taps = root_raised_cosine_pulse(np.arange(-span * sps, span * sps + 1) / sps, 1.0, 0.35) / np.sqrt(sps)
    symbols = bpsk(rng.integers(0, 2, 400))
    upsampled = np.zeros(symbols.size * sps)
    upsampled[::sps] = symbols
    received = np.convolve(np.convolve(upsampled, taps), taps) + rng.normal(0, 0.05, upsampled.size + 2 * taps.size - 2)
    plot_eye_diagram(received, sps, ax=ax2, offset=2 * span * sps + sps // 2 + sps * 10)
    ax2.set_title(r"Eye diagram, root-raised-cosine $\alpha$ = 0.35")

    ebn0_db = np.linspace(0, 20, 200)
    ebn0 = 10 ** (ebn0_db / 10)
    plot_ber_curve(ebn0_db, ber_bpsk_awgn(ebn0), ax=ax3, label="BPSK, AWGN")
    plot_ber_curve(ebn0_db, ser_mqam_awgn(4 * ebn0, 16) / 4, ax=ax3, label="16-QAM, AWGN (approx.)")
    plot_ber_curve(ebn0_db, ber_bpsk_rayleigh(ebn0), ax=ax3, label="BPSK, Rayleigh fading")
    ax3.set(ylim=(1e-7, 1), title="Bit error rates")
    ax3.legend(fontsize=8)


def materials(axes):
    from engineeringkit.materials import crack_growth_life, ramberg_osgood_strain
    from engineeringkit.materials.visualizers import plot_crack_growth, plot_stress_strain, plot_yield_surfaces

    ax1, ax2, ax3 = axes
    for n, label in ((5, "n = 5"), (10, "n = 10"), (25, "n = 25")):
        stress = np.linspace(0, 520e6, 300)
        plot_stress_strain(ramberg_osgood_strain(stress, 200e9, 350e6, n), stress, ax=ax1, label=label)
    ax1.set(xlim=(0, 0.02), title="Ramberg-Osgood stress-strain curves")
    ax1.legend(fontsize=8)

    plot_yield_surfaces(250e6, ax=ax2)
    ax2.set_title("Tresca and von Mises yield surfaces")

    for delta_sigma in (80.0, 100.0, 120.0):
        life = crack_growth_life(1e-11, 3.0, delta_sigma, 1e-3, 25e-3, Y=1.12)
        plot_crack_growth(life.cycles, life.crack_length, ax=ax3, label=f"Δσ = {delta_sigma:g} MPa")
    ax3.set_title("Fatigue crack growth by Paris's law")
    ax3.legend(fontsize=8)


def reliability(axes):
    from engineeringkit.reliability import bathtub_hazard, kaplan_meier, weibull_rank_regression
    from engineeringkit.reliability.visualizers import plot_hazard, plot_survival

    ax1, ax2, ax3 = axes
    ax1.figure.get_layout_engine().set(wspace=0.06)
    t = np.linspace(0.01, 1000, 2000)
    plot_hazard(
        t, bathtub_hazard(t, early_scale=0.02, early_rate=0.03, random_rate=0.002, wear_scale=2e-14, wear_power=4), ax=ax1, log=True
    )
    ax1.set_title("The bathtub curve")

    rng = np.random.default_rng(7)
    times = 1000.0 * rng.weibull(2.2, 30)
    fit = weibull_rank_regression(times)
    fit.plot(ax=ax2)
    ax2.set_title(f"Weibull probability plot: β = {fit.beta:.2f}, η = {fit.eta:.0f} h")

    failure = 1000.0 * rng.weibull(1.5, 60)
    censor = rng.uniform(200, 1500, 60)
    km = kaplan_meier(np.minimum(failure, censor), failure <= censor)
    plot_survival(km.times, km.survival, ax=ax3, label="Kaplan-Meier")
    grid = np.linspace(0, 1500, 300)
    ax3.plot(grid, np.exp(-((grid / 1000.0) ** 1.5)), "k--", label="true Weibull")
    ax3.set_title("Survival from censored life data")
    ax3.legend(fontsize=8)


def hero(axes):
    """The docs landing-page banner: airflow, chaos and a limit cycle."""
    from matplotlib.collections import LineCollection
    from matplotlib.path import Path as MplPath

    from engineeringkit.aerodynamics import naca4_airfoil, panel_method
    from engineeringkit.aerodynamics.systems.panels import _unit_panel_velocities
    from engineeringkit.control_systems import phase_portrait
    from engineeringkit.mechanical_engineering import simulate_lagrangian

    ax1, ax2, ax3 = axes
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Streamlines round a NACA 2412 at 8 degrees: the flow field is rebuilt from
    # the panel method's own source and vortex strengths.
    alpha = np.radians(8.0)
    flow = panel_method(naca4_airfoil("2412", 121), alpha)
    nodes = flow.nodes
    start, delta = nodes[:-1], np.diff(nodes, axis=0)
    length = np.hypot(delta[:, 0], delta[:, 1])
    tangent = delta / length[:, None]
    normal = np.column_stack((-tangent[:, 1], tangent[:, 0]))
    gx, gy = np.meshgrid(np.linspace(-0.7, 1.7, 300), np.linspace(-0.9, 0.9, 225))
    points = np.column_stack((gx.ravel(), gy.ravel()))
    velocity = np.empty_like(points)
    for chunk in np.array_split(np.arange(len(points)), 40):
        source, vortex = _unit_panel_velocities(points[chunk], start, tangent, normal, length)
        velocity[chunk] = [np.cos(alpha), np.sin(alpha)] + np.einsum("ijk,j->ik", source, flow.source_strengths)
        velocity[chunk] += flow.vortex_strength * vortex.sum(axis=1)
    inside = MplPath(nodes).contains_points(points, radius=-1e-3) | MplPath(nodes).contains_points(points, radius=1e-3)
    u = np.ma.masked_array(velocity[:, 0], inside).reshape(gx.shape)
    v = np.ma.masked_array(velocity[:, 1], inside).reshape(gx.shape)
    speed = np.hypot(u, v)
    ax1.contourf(gx, gy, speed, levels=np.linspace(0.0, 1.9, 39), cmap="RdYlBu_r", extend="max")
    seeds = np.column_stack((np.full(36, -0.7), np.linspace(-0.88, 0.88, 36)))
    ax1.streamplot(gx, gy, u, v, color="#1f2430", linewidth=0.6, start_points=seeds, arrowsize=0.6, broken_streamlines=False)
    ax1.fill(nodes[:, 0], nodes[:, 1], color="#3a3f4b", zorder=3)
    ax1.set_aspect("equal")
    ax1.set_xlim(-0.7, 1.7)
    ax1.set_ylim(-0.9, 0.9)
    ax1.set_title(f"Panel-method flow over a NACA 2412,  $C_L$ = {flow.cl:.2f}")

    # Two double pendulums released 0.001 rad apart part company within seconds.
    g = 9.81

    def lagrangian(q, qd):
        kinetic = qd[0] ** 2 + 0.5 * qd[1] ** 2 + qd[0] * qd[1] * np.cos(q[0] - q[1])
        return kinetic + g * (2 * np.cos(q[0]) + np.cos(q[1]))

    t = np.linspace(0.0, 18.0, 5400)
    for offset, cmap in ((0.0, "Blues"), (1e-3, "Oranges")):
        q = simulate_lagrangian(lagrangian, [2.4, 2.4 + offset], [0.0, 0.0], t).displacement
        tip = np.column_stack((np.sin(q[:, 0]) + np.sin(q[:, 1]), -np.cos(q[:, 0]) - np.cos(q[:, 1])))
        segments = np.stack((tip[:-1], tip[1:]), axis=1)
        trail = LineCollection(segments, cmap=cmap, linewidth=0.8, alpha=0.75)
        trail.set_array(t[:-1])
        trail.set_clim(-0.3 * t[-1], t[-1])
        ax2.add_collection(trail)
        elbow = [np.sin(q[-1, 0]), -np.cos(q[-1, 0])]
        ax2.plot([0, elbow[0], tip[-1, 0]], [0, elbow[1], tip[-1, 1]], "-o", color=plt.get_cmap(cmap)(1.0), lw=2, ms=6)
    ax2.plot(0, 0, "o", color="black", ms=5)
    ax2.set_aspect("equal")
    ax2.set_xlim(-2.2, 2.2)
    ax2.set_ylim(-2.2, 2.2)
    ax2.set_title("Lagrangian double pendulum: 0.001 rad apart, chaos")

    # The Van der Pol oscillator: every start winds onto the same limit cycle.
    mu = 1.5

    def van_der_pol(x, v):
        return mu * (1 - x**2) * v - x

    gx, gv = np.meshgrid(np.linspace(-4, 4, 200), np.linspace(-6, 6, 200))
    field = van_der_pol(gx, gv)
    ax3.streamplot(gx, gv, gv, field, color=np.log1p(np.hypot(gv, field)), cmap="viridis", linewidth=0.7, density=1.6, arrowsize=0.7)
    t = np.linspace(0.0, 30.0, 6000)
    starts = [(0.1, 0.0), (-3.0, 2.0), (3.0, -2.0), (-1.0, -4.5), (1.0, 4.5)]
    for x0, v0 in starts:
        path = phase_portrait(van_der_pol, x0, v0, t)
        ax3.plot(path.x, path.v, color="#f77f00", lw=1.4, alpha=0.9)
        ax3.plot(x0, v0, "o", color="#f77f00", ms=6)
    cycle = phase_portrait(van_der_pol, 2.0, 0.0, np.linspace(0.0, 40.0, 8000))
    ax3.plot(cycle.x[-1200:], cycle.v[-1200:], color="#c1121f", lw=2.6)
    ax3.set_xlim(-4, 4)
    ax3.set_ylim(-6, 6)
    ax3.set_title(rf"Van der Pol limit cycle, $\mu$ = {mu}")


FIGURES = {
    f.__name__: f
    for f in (
        structural_engineering,
        mechanical_engineering,
        thermodynamics,
        fluid_engineering,
        heat_transfer,
        aerodynamics,
        electrical_engineering,
        control_systems,
        signal_processing,
        communications,
        materials,
        reliability,
        hero,
    )
}


def main(names):
    for name in names or FIGURES:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
        FIGURES[name](axes)
        out = OUT_DIR / f"readme_{name}.png"
        fig.savefig(out, dpi=110)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
