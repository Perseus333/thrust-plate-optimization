import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

BORE_RADIUS = 0.029
LOAD_OUTER_RADIUS = 0.036
CLAMP_RADIUS = 0.049
PLATE_OUTER_RADIUS = 0.055

THICKNESS = 0.006
OUTER_DIAMETER = 0.110

YOUNGS_MODULUS = 70.0e9
POISSON_RATIO = 0.33
DENSITY = 2700.0
YIELD_STRESS = 250.0e6

PEAK_THRUST = 1261.0
SAFETY_FACTOR_REQUIRED = 2.5
ULTIMATE_THRUST = PEAK_THRUST * SAFETY_FACTOR_REQUIRED

ALLOWABLE_STRESS = YIELD_STRESS / SAFETY_FACTOR_REQUIRED

SAMPLE_POINTS = 401

PLATE_AREA = math.pi * (PLATE_OUTER_RADIUS**2 - BORE_RADIUS**2)
LOAD_AREA = math.pi * (LOAD_OUTER_RADIUS**2 - BORE_RADIUS**2)

BASIS = {
    "deflection": lambda r: np.array([1.0, r**2, math.log(r), r**2 * math.log(r)]),
    "slope": lambda r: np.array([0.0, 2.0 * r, 1.0 / r, 2.0 * r * math.log(r) + r]),
    "curvature": lambda r: np.array([0.0, 2.0, -1.0 / r**2, 2.0 * math.log(r) + 3.0]),
    "shear": lambda r: np.array([0.0, 0.0, 0.0, 1.0]),
}

TERMS = {
    "deflection": lambda r: r**4 / 64.0,
    "slope": lambda r: r**3 / 16.0,
    "curvature": lambda r: 3.0 * r**2 / 16.0,
    "radial_moment": lambda r: -(3.0 + POISSON_RATIO) * r**2 / 16.0,
    "hoop_moment": lambda r: -(1.0 + 3.0 * POISSON_RATIO) * r**2 / 16.0,
    "shear": lambda r: -r / 2.0,
}


def stiffness(thickness):
    return YOUNGS_MODULUS * thickness**3 / (12.0 * (1.0 - POISSON_RATIO**2))


def homogeneous(r, rigidity, kind):
    if kind == "radial_moment":
        return -rigidity * (
            BASIS["curvature"](r) + POISSON_RATIO * BASIS["slope"](r) / r
        )
    if kind == "hoop_moment":
        return -rigidity * (
            BASIS["slope"](r) / r + POISSON_RATIO * BASIS["curvature"](r)
        )
    if kind == "shear":
        return -4.0 * rigidity * BASIS["shear"](r) / r
    return BASIS[kind](r)


def particular(region, r, pressure, rigidity, kind):
    if region != 0:
        return 0.0
    if kind == "deflection":
        return pressure * TERMS["deflection"](r) / rigidity
    if kind == "slope":
        return pressure * TERMS["slope"](r) / rigidity
    if kind == "curvature":
        return pressure * TERMS["curvature"](r) / rigidity
    return pressure * TERMS[kind](r)


def coefficients(thickness=THICKNESS, clamp_radius=CLAMP_RADIUS, total_force=PEAK_THRUST):
    pressure = total_force / LOAD_AREA
    rigidity = stiffness(thickness)
    joint = LOAD_OUTER_RADIUS
    rows = [
        (0, joint, "deflection", True),
        (0, joint, "slope", True),
        (0, joint, "radial_moment", True),
        (0, joint, "shear", True),
        (0, BORE_RADIUS, "radial_moment", False),
        (0, BORE_RADIUS, "shear", False),
        (1, clamp_radius, "deflection", False),
        (1, clamp_radius, "slope", False),
    ]
    matrix = np.zeros((len(rows), len(rows)))
    right = np.zeros(len(rows))
    for row, (region, r, kind, shared) in enumerate(rows):
        matrix[row, region * 4 : region * 4 + 4] = homogeneous(r, rigidity, kind)
        if shared:
            matrix[row, 4:8] = -homogeneous(r, rigidity, kind)
        right[row] = -particular(region, r, pressure, rigidity, kind)
        if shared:
            right[row] += particular(1, r, pressure, rigidity, kind)
    return np.linalg.solve(matrix, right)


def solve(thickness=THICKNESS, clamp_radius=CLAMP_RADIUS, total_force=PEAK_THRUST):
    if not BORE_RADIUS < clamp_radius <= PLATE_OUTER_RADIUS:
        raise ValueError("clamp radius must lie inside the plate")
    return {
        "coefficients": coefficients(thickness, clamp_radius, total_force),
        "rigidity": stiffness(thickness),
        "thickness": thickness,
        "pressure": total_force / LOAD_AREA,
        "clamp_radius": clamp_radius,
    }


def weight(solution, r):
    if r < LOAD_OUTER_RADIUS:
        return 1.0
    return 0.0


def term_at(solution, r, kind):
    rigidity = solution["rigidity"]
    pressure = solution["pressure"] * weight(solution, r)
    region = 0 if r < LOAD_OUTER_RADIUS else 1
    local = solution["coefficients"][region * 4 : region * 4 + 4]
    if kind in ("deflection", "slope", "curvature"):
        return float(homogeneous(r, rigidity, kind) @ local) + (
            pressure * TERMS[kind](r) / rigidity
        )
    return float(homogeneous(r, rigidity, kind) @ local) + pressure * TERMS[kind](r)


def field(solution, r, kind):
    rigidity = solution["rigidity"]
    if kind == "radial_stress":
        return 6.0 * term_at(solution, r, "radial_moment") / solution["thickness"] ** 2
    if kind == "hoop_stress":
        return 6.0 * term_at(solution, r, "hoop_moment") / solution["thickness"] ** 2
    if kind == "von_mises":
        radial = field(solution, r, "radial_stress")
        hoop = field(solution, r, "hoop_stress")
        return math.sqrt(radial**2 - radial * hoop + hoop**2)
    if kind == "shear":
        step = 1.0e-9
        ahead = term_at(solution, r + step, "curvature") + term_at(
            solution, r + step, "slope"
        ) / (r + step)
        behind = term_at(solution, r - step, "curvature") + term_at(
            solution, r - step, "slope"
        ) / (r - step)
        return -rigidity * (ahead - behind) / (2.0 * step)
    return term_at(solution, r, kind)


def peak(solution, kind, r_start=BORE_RADIUS, r_end=CLAMP_RADIUS):
    radii = np.linspace(r_start, r_end, SAMPLE_POINTS)
    values = np.array([field(solution, r, kind) for r in radii])
    index = int(np.argmax(np.abs(values)))
    return float(radii[index]), float(values[index])


def profile(solution, r_start=BORE_RADIUS, r_end=CLAMP_RADIUS):
    radii = np.linspace(r_start, r_end, SAMPLE_POINTS)
    return radii, {
        kind: np.array([field(solution, r, kind) for r in radii])
        for kind in ("radial_stress", "hoop_stress", "von_mises", "deflection")
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
        solution = solve(total_force=force)
        radius, stress = peak(solution, "von_mises")
        _, deflection = peak(solution, "deflection")
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
    print(f"Optimised for safety factor {SAFETY_FACTOR_REQUIRED:.2f} at {ULTIMATE_THRUST:.0f} N")
    print("  thickness      {:8.4f} mm".format(result["thickness"] * 1e3))
    print("  peak stress    {:8.4f} MPa".format(result["peak_von_mises_stress"] / 1e6))
    print("  safety factor  {:8.4f}".format(result["safety_factor"]))
    print("  deflection     {:8.4f} um".format(result["max_deflection"] * 1e6))
    print("  mass           {:8.4f} g".format(result["mass"] * 1e3))
    print("  mass saving    {:8.4f} %".format(result["mass_saving_fraction"] * 100.0))


def optimise(total_force=ULTIMATE_THRUST, clamp_radius=CLAMP_RADIUS):
    reference = peak(solve(clamp_radius=clamp_radius, total_force=total_force), "von_mises")[1]
    thickness = THICKNESS * math.sqrt(reference / ALLOWABLE_STRESS)
    solution = solve(thickness=thickness, clamp_radius=clamp_radius, total_force=total_force)
    radius, stress = peak(solution, "von_mises")
    _, deflection = peak(solution, "deflection")
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


def plot(path="stress_profile.png", clamp_radius=CLAMP_RADIUS):
    as_built = solve(clamp_radius=clamp_radius)
    radii, curves = profile(as_built, r_end=clamp_radius)
    thickness = THICKNESS * math.sqrt(
        field(as_built, clamp_radius, "von_mises") / ALLOWABLE_STRESS
    )
    optimised = solve(thickness=thickness, clamp_radius=clamp_radius, total_force=ULTIMATE_THRUST)
    figure, axes = plt.subplots(3, 1, figsize=(7.0, 11.0))
    axes[0].plot(radii * 1e3, curves["radial_stress"] / 1e6, label="radial stress")
    axes[0].plot(radii * 1e3, curves["hoop_stress"] / 1e6, label="hoop stress")
    axes[0].plot(radii * 1e3, curves["von_mises"] / 1e6, label="von Mises", linestyle="--")
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
        radii * 1e3, curves["von_mises"] / 1e6, label=f"{THICKNESS * 1e3:.2f} mm as built"
    )
    axes[2].plot(
        radii * 1e3,
        np.array([field(optimised, r, "von_mises") for r in radii]) / 1e6,
        label=f"{optimised['thickness'] * 1e3:.2f} mm optimised",
    )
    axes[2].axhline(YIELD_STRESS / 1e6, color="red", linestyle=":", label="yield 6082-T6")
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
