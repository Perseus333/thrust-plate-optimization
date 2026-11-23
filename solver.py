import math

import numpy as np

from model import (
    BORE_RADIUS,
    CLAMP_RADIUS,
    LOAD_AREA,
    LOAD_OUTER_RADIUS,
    PEAK_THRUST,
    PLATE_OUTER_RADIUS,
    POISSON_RATIO,
    SAMPLE_POINTS,
    THICKNESS,
    YOUNGS_MODULUS,
)

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
    if kind == "radial_stress":
        return 6.0 * term_at(solution, r, "radial_moment") / solution["thickness"] ** 2
    if kind == "hoop_stress":
        return 6.0 * term_at(solution, r, "hoop_moment") / solution["thickness"] ** 2
    if kind == "von_mises":
        radial = field(solution, r, "radial_stress")
        hoop = field(solution, r, "hoop_stress")
        return math.sqrt(radial**2 - radial * hoop + hoop**2)
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
