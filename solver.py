import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


BORE_RADIUS = 0.029
LOAD_INNER_RADIUS = 0.029
LOAD_OUTER_RADIUS = 0.036
FREE_OUTER_RADIUS = 0.043
PLATE_OUTER_RADIUS = 0.055
PLATE_INNER_RADIUS = 0.029

THICKNESS = 0.006
OUTER_DIAMETER = 0.110

YOUNGS_MODULUS = 70.0e9
POISSON_RATIO = 0.33
DENSITY = 2700.0
YIELD_STRESS = 250.0e6
ULTIMATE_STRESS = 300.0e6

PEAK_THRUST = 1261.0
SAFETY_FACTOR_REQUIRED = 2.5
ULTIMATE_THRUST = PEAK_THRUST * SAFETY_FACTOR_REQUIRED
BURN_TIME = 2.7
TOTAL_IMPULSE = 2714.0

ALLOWABLE_STRESS = YIELD_STRESS / SAFETY_FACTOR_REQUIRED

BASIS_COUNT = 4
SAMPLE_POINTS = 401


def flexural_rigidity(thickness):
    return YOUNGS_MODULUS * thickness**3 / (12.0 * (1.0 - POISSON_RATIO**2))


def annulus_area(r_inner, r_outer):
    return math.pi * (r_outer**2 - r_inner**2)


def plate_volume(thickness):
    return annulus_area(PLATE_INNER_RADIUS, PLATE_OUTER_RADIUS) * thickness


def plate_mass(thickness):
    return DENSITY * plate_volume(thickness)


def load_area():
    return annulus_area(LOAD_INNER_RADIUS, LOAD_OUTER_RADIUS)


def uniform_pressure(total_force):
    return total_force / load_area()


def load_bands(total_force):
    return [(LOAD_INNER_RADIUS, LOAD_OUTER_RADIUS, uniform_pressure(total_force))]


def basis(r):
    return np.array([1.0, r**2, math.log(r), r**2 * math.log(r)])


def basis_slope(r):
    return np.array([0.0, 2.0 * r, 1.0 / r, 2.0 * r * math.log(r) + r])


def basis_curvature(r):
    return np.array([0.0, 2.0, -1.0 / r**2, 2.0 * math.log(r) + 3.0])


def basis_shear(r, rigidity):
    return (-rigidity) * np.array([0.0, 0.0, 0.0, 4.0 / r])


def basis_radial_moment(r, rigidity, poisson_ratio):
    return (-rigidity) * (basis_curvature(r) + poisson_ratio * basis_slope(r) / r)


def basis_hoop_moment(r, rigidity, poisson_ratio):
    return (-rigidity) * (basis_slope(r) / r + poisson_ratio * basis_curvature(r))


def particular_deflection(r, pressure, rigidity):
    return pressure * r**4 / (64.0 * rigidity)


def particular_slope(r, pressure, rigidity):
    return pressure * r**3 / (16.0 * rigidity)


def particular_curvature(r, pressure, rigidity):
    return 3.0 * pressure * r**2 / (16.0 * rigidity)


def particular_shear(r, pressure):
    return -pressure * r / 2.0


def particular_radial_moment(r, pressure, poisson_ratio):
    return -(3.0 + poisson_ratio) * pressure * r**2 / 16.0


def particular_hoop_moment(r, pressure, poisson_ratio):
    return -(1.0 + 3.0 * poisson_ratio) * pressure * r**2 / 16.0


def pressure_at(bands, r):
    for r_inner, r_outer, pressure in bands:
        if r_inner <= r <= r_outer:
            return pressure
    return 0.0


def collect_cuts(bands, r_inner, r_outer):
    cuts = {r_inner, r_outer}
    for band_inner, band_outer, _ in bands:
        if r_inner < band_inner < r_outer:
            cuts.add(band_inner)
        if r_inner < band_outer < r_outer:
            cuts.add(band_outer)
    return sorted(cuts)


def build_regions(bands, r_inner, r_outer):
    edges = collect_cuts(bands, r_inner, r_outer)
    regions = []
    for index in range(len(edges) - 1):
        r_start = edges[index]
        r_end = edges[index + 1]
        midpoint = 0.5 * (r_start + r_end)
        regions.append((r_start, r_end, pressure_at(bands, midpoint)))
    return regions


def region_index(regions, r):
    if r <= regions[0][0]:
        return 0
    if r >= regions[-1][1]:
        return len(regions) - 1
    for index, (r_start, r_end, _) in enumerate(regions):
        if r_start <= r < r_end:
            return index
    raise ValueError("radius lies outside the solved domain")


class Edge:
    def __init__(self, radius, condition):
        self.radius = radius
        self.condition = condition


class PlateSolution:
    def __init__(self, regions, coefficients, rigidity, poisson_ratio, thickness):
        self.regions = regions
        self.coefficients = coefficients
        self.rigidity = rigidity
        self.poisson_ratio = poisson_ratio
        self.thickness = thickness

    def local(self, r):
        index = region_index(self.regions, r)
        start = index * BASIS_COUNT
        return self.regions[index][2], self.coefficients[start:start + BASIS_COUNT]

    def deflection(self, r):
        pressure, local = self.local(r)
        return float(basis(r) @ local) + particular_deflection(r, pressure, self.rigidity)

    def slope(self, r):
        pressure, local = self.local(r)
        return float(basis_slope(r) @ local) + particular_slope(r, pressure, self.rigidity)

    def curvature(self, r):
        pressure, local = self.local(r)
        return float(basis_curvature(r) @ local) + particular_curvature(r, pressure, self.rigidity)

    def radial_moment(self, r):
        pressure, local = self.local(r)
        homogeneous = float(basis_radial_moment(r, self.rigidity, self.poisson_ratio) @ local)
        return homogeneous + particular_radial_moment(r, pressure, self.poisson_ratio)

    def hoop_moment(self, r):
        pressure, local = self.local(r)
        homogeneous = float(basis_hoop_moment(r, self.rigidity, self.poisson_ratio) @ local)
        return homogeneous + particular_hoop_moment(r, pressure, self.poisson_ratio)

    def radial_stress(self, r):
        return 6.0 * self.radial_moment(r) / self.thickness**2

    def hoop_stress(self, r):
        return 6.0 * self.hoop_moment(r) / self.thickness**2

    def von_mises_stress(self, r):
        radial = self.radial_stress(r)
        hoop = self.hoop_stress(r)
        return math.sqrt(radial**2 - radial * hoop + hoop**2)

    def shear(self, r, step=1.0e-9):
        combined = self.curvature(r) + self.slope(r) / r
        ahead = self.curvature(r + step) + self.slope(r + step) / (r + step)
        behind = self.curvature(r - step) + self.slope(r - step) / (r - step)
        return -self.rigidity * (ahead - behind) / (2.0 * step)


def assemble(regions, inner_edge, outer_edge, rigidity, poisson_ratio):
    count = BASIS_COUNT * len(regions)
    matrix = np.zeros((count, count))
    right = np.zeros(count)
    row = 0

    def slice_of(region):
        return slice(region * BASIS_COUNT, (region + 1) * BASIS_COUNT)

    def add(region, sign, vector):
        matrix[row, slice_of(region)] += sign * vector

    def value(region, r, kind):
        pressure = regions[region][2]
        if kind == "w":
            return particular_deflection(r, pressure, rigidity)
        if kind == "slope":
            return particular_slope(r, pressure, rigidity)
        if kind == "shear":
            return particular_shear(r, pressure)
        if kind == "radial_moment":
            return particular_radial_moment(r, pressure, poisson_ratio)
        raise ValueError("unknown particular solution")

    def add_edge_condition(region, r, condition):
        nonlocal row
        if condition == "clamped":
            matrix[row] = 0.0
            add(region, 1.0, basis(r))
            right[row] = -value(region, r, "w")
            row += 1
            matrix[row] = 0.0
            add(region, 1.0, basis_slope(r))
            right[row] = -value(region, r, "slope")
            row += 1
        elif condition == "free":
            matrix[row] = 0.0
            add(region, 1.0, basis_radial_moment(r, rigidity, poisson_ratio))
            right[row] = -value(region, r, "radial_moment")
            row += 1
            matrix[row] = 0.0
            add(region, 1.0, basis_shear(r, rigidity))
            right[row] = -value(region, r, "shear")
            row += 1
        else:
            raise ValueError("unknown edge condition")

    for region in range(len(regions) - 1):
        shared = regions[region][1]
        other = region + 1

        matrix[row] = 0.0
        add(region, 1.0, basis(shared))
        add(other, -1.0, basis(shared))
        right[row] = -(value(region, shared, "w") - value(other, shared, "w"))
        row += 1

        matrix[row] = 0.0
        add(region, 1.0, basis_slope(shared))
        add(other, -1.0, basis_slope(shared))
        right[row] = -(value(region, shared, "slope") - value(other, shared, "slope"))
        row += 1

        matrix[row] = 0.0
        add(region, 1.0, basis_radial_moment(shared, rigidity, poisson_ratio))
        add(other, -1.0, basis_radial_moment(shared, rigidity, poisson_ratio))
        right[row] = -(value(region, shared, "radial_moment") - value(other, shared, "radial_moment"))
        row += 1

        matrix[row] = 0.0
        add(region, 1.0, basis_shear(shared, rigidity))
        add(other, -1.0, basis_shear(shared, rigidity))
        right[row] = -(value(region, shared, "shear") - value(other, shared, "shear"))
        row += 1

    add_edge_condition(0, regions[0][0], inner_edge.condition)
    add_edge_condition(len(regions) - 1, regions[-1][1], outer_edge.condition)

    return matrix, right


def solve(loads, inner_edge, outer_edge, thickness=None):
    if thickness is None:
        thickness = THICKNESS
    if inner_edge.radius >= outer_edge.radius:
        raise ValueError("inner edge radius must be smaller than outer edge radius")
    for band_inner, band_outer, _ in loads:
        if band_inner < inner_edge.radius or band_outer > outer_edge.radius:
            raise ValueError("load band lies outside the solved domain")
    rigidity = flexural_rigidity(thickness)
    regions = build_regions(loads, inner_edge.radius, outer_edge.radius)
    matrix, right = assemble(regions, inner_edge, outer_edge, rigidity, POISSON_RATIO)
    coefficients = np.linalg.solve(matrix, right)
    return PlateSolution(regions, coefficients, rigidity, POISSON_RATIO, thickness)


def solve_layout(total_force, clamp_radius, thickness=None):
    inner = Edge(BORE_RADIUS, "free")
    outer = Edge(clamp_radius, "clamped")
    return solve(load_bands(total_force), inner, outer, thickness=thickness)


def sample(solution, r_start, r_end, points=SAMPLE_POINTS):
    radii = np.linspace(r_start, r_end, points)
    radial_stress = np.array([solution.radial_stress(r) for r in radii])
    hoop_stress = np.array([solution.hoop_stress(r) for r in radii])
    von_mises = np.array([solution.von_mises_stress(r) for r in radii])
    deflection = np.array([solution.deflection(r) for r in radii])
    return radii, radial_stress, hoop_stress, von_mises, deflection


def summary(solution, r_start, r_end, points=SAMPLE_POINTS):
    radii, radial_stress, hoop_stress, von_mises, deflection = sample(
        solution, r_start, r_end, points
    )
    peak = int(np.argmax(von_mises))
    radial_peak = int(np.argmax(np.abs(radial_stress)))
    hoop_peak = int(np.argmax(np.abs(hoop_stress)))
    deflection_peak = int(np.argmax(np.abs(deflection)))
    return {
        "peak_von_mises_stress": float(von_mises[peak]),
        "peak_von_mises_radius": float(radii[peak]),
        "peak_radial_stress": float(radial_stress[radial_peak]),
        "peak_hoop_stress": float(hoop_stress[hoop_peak]),
        "max_deflection": float(deflection[deflection_peak]),
        "max_deflection_radius": float(radii[deflection_peak]),
        "safety_factor": float(YIELD_STRESS / von_mises[peak]),
    }


def required_thickness(reference_stress, allowable_stress, reference_thickness=THICKNESS):
    return reference_thickness * math.sqrt(reference_stress / allowable_stress)


def optimise_for_safety_factor(total_force=ULTIMATE_THRUST, clamp_radius=FREE_OUTER_RADIUS):
    reference = solve_layout(total_force, clamp_radius, thickness=THICKNESS)
    reference_stress = summary(reference, BORE_RADIUS, clamp_radius)["peak_von_mises_stress"]
    thickness = required_thickness(reference_stress, ALLOWABLE_STRESS)
    optimised = solve_layout(total_force, clamp_radius, thickness=thickness)
    result = summary(optimised, BORE_RADIUS, clamp_radius)
    result["thickness"] = thickness
    result["mass"] = plate_mass(thickness)
    result["baseline_mass"] = plate_mass(THICKNESS)
    result["mass_saving"] = result["baseline_mass"] - result["mass"]
    result["mass_saving_fraction"] = result["mass_saving"] / result["baseline_mass"]
    return result


def print_summary(clamp_radius=FREE_OUTER_RADIUS):
    print("Annular thrust plate, Al 6082-T6")
    print("Bore {:.0f} mm, load band {:.0f}-{:.0f} mm, clamped {:.0f}-{:.0f} mm".format(
        BORE_RADIUS * 1e3,
        LOAD_INNER_RADIUS * 1e3,
        LOAD_OUTER_RADIUS * 1e3,
        clamp_radius * 1e3,
        PLATE_OUTER_RADIUS * 1e3,
    ))
    print("Load area {:.1f} mm^2, pressure at {:.0f} N {:.1f} kPa".format(
        load_area() * 1e6, PEAK_THRUST, uniform_pressure(PEAK_THRUST) / 1e3
    ))
    print()
    for force, label in [(PEAK_THRUST, "peak"), (ULTIMATE_THRUST, "ultimate")]:
        solution = solve_layout(force, clamp_radius)
        result = summary(solution, BORE_RADIUS, clamp_radius)
        print("{:9s} {:7.0f} N  peak {:8.4f} MPa at r={:5.1f} mm  {:9.4f} um  SF {:6.2f}".format(
            label,
            force,
            result["peak_von_mises_stress"] / 1e6,
            result["peak_von_mises_radius"] * 1e3,
            result["max_deflection"] * 1e6,
            result["safety_factor"],
        ))
    print()
    result = optimise_for_safety_factor()
    print("Optimised for safety factor {:.2f} at {:.0f} N".format(
        SAFETY_FACTOR_REQUIRED, ULTIMATE_THRUST
    ))
    print("  thickness      {:8.4f} mm".format(result["thickness"] * 1e3))
    print("  peak stress    {:8.4f} MPa".format(result["peak_von_mises_stress"] / 1e6))
    print("  safety factor  {:8.4f}".format(result["safety_factor"]))
    print("  deflection     {:8.4f} um".format(result["max_deflection"] * 1e6))
    print("  mass           {:8.4f} g".format(result["mass"] * 1e3))
    print("  mass saving    {:8.4f} %".format(result["mass_saving_fraction"] * 100.0))


def plot(clamp_radius=FREE_OUTER_RADIUS, path="stress_profile.png"):
    solution = solve_layout(PEAK_THRUST, clamp_radius)
    radii, radial_stress, hoop_stress, von_mises, deflection = sample(
        solution, BORE_RADIUS, clamp_radius
    )
    optimised = optimise_for_safety_factor(clamp_radius=clamp_radius)
    optimised_solution = solve_layout(ULTIMATE_THRUST, clamp_radius, thickness=optimised["thickness"])
    optimised_radii, _, _, optimised_von_mises, _ = sample(
        optimised_solution, BORE_RADIUS, clamp_radius
    )

    millimetres = radii * 1e3
    figure, axes = plt.subplots(3, 1, figsize=(7.0, 11.0))

    axes[0].plot(millimetres, radial_stress / 1e6, label="radial stress")
    axes[0].plot(millimetres, hoop_stress / 1e6, label="hoop stress")
    axes[0].plot(millimetres, von_mises / 1e6, label="von Mises", linestyle="--")
    axes[0].axhline(YIELD_STRESS / 1e6, color="red", linestyle=":", label="yield 6082-T6")
    axes[0].axvline(LOAD_OUTER_RADIUS * 1e3, color="grey", linestyle="-.", label="load band edge")
    axes[0].set_xlabel("radius, mm")
    axes[0].set_ylabel("bending stress, MPa")
    axes[0].set_title("Stress at {:.0f} N, t = {:.2f} mm".format(PEAK_THRUST, THICKNESS * 1e3))
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8)

    axes[1].plot(millimetres, deflection * 1e6, color="black")
    axes[1].set_xlabel("radius, mm")
    axes[1].set_ylabel("deflection, um")
    axes[1].set_title("Transverse deflection at {:.0f} N".format(PEAK_THRUST))
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(radii * 1e3, von_mises / 1e6, label="{:.2f} mm as built".format(THICKNESS * 1e3))
    axes[2].plot(
        optimised_radii * 1e3,
        optimised_von_mises / 1e6,
        label="{:.2f} mm optimised".format(optimised["thickness"] * 1e3),
    )
    axes[2].axhline(YIELD_STRESS / 1e6, color="red", linestyle=":", label="yield 6082-T6")
    axes[2].set_xlabel("radius, mm")
    axes[2].set_ylabel("von Mises stress, MPa")
    axes[2].set_title("von Mises at the {:.0f} N ultimate load".format(ULTIMATE_THRUST))
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def main():
    print_summary()
    print()
    print("wrote " + plot())


if __name__ == "__main__":
    main()
