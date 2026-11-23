import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import solver
from model import (
    ALLOWABLE_STRESS,
    BORE_RADIUS,
    CLAMP_RADIUS,
    DENSITY,
    LOAD_AREA,
    LOAD_OUTER_RADIUS,
    PEAK_THRUST,
    PLATE_AREA,
    PLATE_OUTER_RADIUS,
    SAFETY_FACTOR_REQUIRED,
    THICKNESS,
    ULTIMATE_THRUST,
    YIELD_STRESS,
)


def optimise(total_force=ULTIMATE_THRUST, clamp_radius=CLAMP_RADIUS):
    reference = solver.peak(
        solver.solve(clamp_radius=clamp_radius, total_force=total_force), "von_mises"
    )[1]
    thickness = THICKNESS * math.sqrt(reference / ALLOWABLE_STRESS)
    solution = solver.solve(
        thickness=thickness, clamp_radius=clamp_radius, total_force=total_force
    )
    radius, stress = solver.peak(solution, "von_mises")
    _, deflection = solver.peak(solution, "deflection")
    mass = DENSITY * PLATE_AREA * thickness
    baseline = DENSITY * PLATE_AREA * THICKNESS
    return {
        "thickness": thickness,
        "peak_von_mises_stress": stress,
        "peak_von_mises_radius": radius,
        "max_deflection": deflection,
        "safety_factor": YIELD_STRESS / stress,
        "mass": mass,
        "baseline_mass": baseline,
        "mass_saving": baseline - mass,
        "mass_saving_fraction": (baseline - mass) / baseline,
    }


def report():
    print("Annular thrust plate, Al 6082-T6")
    print(
        f"Bore {BORE_RADIUS * 1e3:.0f} mm, load band {BORE_RADIUS * 1e3:.0f}-{LOAD_OUTER_RADIUS * 1e3:.0f} mm, clamped {CLAMP_RADIUS * 1e3:.0f}-{PLATE_OUTER_RADIUS * 1e3:.0f} mm"
    )
    print(
        f"Load area {LOAD_AREA * 1e6:.1f} mm^2, pressure at {PEAK_THRUST:.0f} N {PEAK_THRUST / LOAD_AREA / 1e3:.1f} kPa"
    )
    print()
    for force, label in [(PEAK_THRUST, "peak"), (ULTIMATE_THRUST, "ultimate")]:
        solution = solver.solve(total_force=force)
        radius, stress = solver.peak(solution, "von_mises")
        _, deflection = solver.peak(solution, "deflection")
        print(
            "{:9s} {:7.0f} N  peak {:8.4f} MPa at r={:5.1f} mm  {:9.4f} um  SF {:6.2f}".format(
                label,
                force,
                stress / 1e6,
                radius * 1e3,
                deflection * 1e6,
                YIELD_STRESS / stress,
            )
        )
    print()
    result = optimise()
    print(
        f"Optimised for safety factor {SAFETY_FACTOR_REQUIRED:.2f} at {ULTIMATE_THRUST:.0f} N"
    )
    print("  thickness      {:8.4f} mm".format(result["thickness"] * 1e3))
    print("  peak stress    {:8.4f} MPa".format(result["peak_von_mises_stress"] / 1e6))
    print("  safety factor  {:8.4f}".format(result["safety_factor"]))
    print("  deflection     {:8.4f} um".format(result["max_deflection"] * 1e6))
    print("  mass           {:8.4f} g".format(result["mass"] * 1e3))
    print("  mass saving    {:8.4f} %".format(result["mass_saving_fraction"] * 100.0))


def plot(path="stress_profile.png", clamp_radius=CLAMP_RADIUS):
    as_built = solver.solve(clamp_radius=clamp_radius)
    radii, curves = solver.profile(as_built, r_end=clamp_radius)
    thickness = THICKNESS * math.sqrt(
        solver.field(as_built, clamp_radius, "von_mises") / ALLOWABLE_STRESS
    )
    optimised = solver.solve(
        thickness=thickness, clamp_radius=clamp_radius, total_force=ULTIMATE_THRUST
    )
    figure, axes = plt.subplots(3, 1, figsize=(7.0, 11.0))
    axes[0].plot(radii * 1e3, curves["radial_stress"] / 1e6, label="radial stress")
    axes[0].plot(radii * 1e3, curves["hoop_stress"] / 1e6, label="hoop stress")
    axes[0].plot(
        radii * 1e3, curves["von_mises"] / 1e6, label="von Mises", linestyle="--"
    )
    axes[0].axhline(YIELD_STRESS / 1e6, color="red", linestyle=":", label="yield 6082-T6")
    axes[0].axvline(
        LOAD_OUTER_RADIUS * 1e3, color="grey", linestyle="-.", label="load band edge"
    )
    axes[0].set_xlabel("radius, mm")
    axes[0].set_ylabel("bending stress, MPa")
    axes[0].set_title(f"Stress at {PEAK_THRUST:.0f} N, t = {THICKNESS * 1e3:.2f} mm")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].plot(radii * 1e3, curves["deflection"] * 1e6, color="black")
    axes[1].set_xlabel("radius, mm")
    axes[1].set_ylabel("deflection, um")
    axes[1].set_title(f"Transverse deflection at {PEAK_THRUST:.0f} N")
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(
        radii * 1e3,
        curves["von_mises"] / 1e6,
        label=f"{THICKNESS * 1e3:.2f} mm as built",
    )
    axes[2].plot(
        radii * 1e3,
        np.array([solver.field(optimised, r, "von_mises") for r in radii]) / 1e6,
        label=f"{optimised['thickness'] * 1e3:.2f} mm optimised",
    )
    axes[2].axhline(
        YIELD_STRESS / 1e6, color="red", linestyle=":", label="yield 6082-T6"
    )
    axes[2].set_xlabel("radius, mm")
    axes[2].set_ylabel("von Mises stress, MPa")
    axes[2].set_title(f"von Mises at the {ULTIMATE_THRUST:.0f} N ultimate load")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def main():
    report()
    print()
    print("wrote " + plot())


if __name__ == "__main__":
    main()
