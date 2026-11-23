import math

import numpy as np

import geometry
import material


BASIS_COUNT = 4


def basis(r):
    return np.array([1.0, r**2, math.log(r), r**2 * math.log(r)])


def basis_slope(r):
    return np.array([0.0, 2.0 * r, 1.0 / r, 2.0 * r * math.log(r) + r])


def basis_curvature(r):
    return np.array([0.0, 2.0, -1.0 / r**2, 2.0 * math.log(r) + 3.0])


def basis_shear(r, rigidity):
    return (-rigidity) * np.array([0.0, 0.0, 0.0, 4.0 / r])


def basis_radial_moment(r, rigidity, poisson_ratio):
    slope = basis_slope(r)
    curvature = basis_curvature(r)
    return (-rigidity) * (curvature + poisson_ratio * slope / r)


def basis_hoop_moment(r, rigidity, poisson_ratio):
    slope = basis_slope(r)
    curvature = basis_curvature(r)
    return (-rigidity) * (slope / r + poisson_ratio * curvature)


class LoadBand:
    def __init__(self, r_inner, r_outer, pressure):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.pressure = pressure


def make_bands(loads):
    bands = []
    for r_inner, r_outer, pressure in loads:
        bands.append(LoadBand(r_inner, r_outer, pressure))
    return bands


def pressure_at(bands, r):
    for band in bands:
        if band.r_inner <= r <= band.r_outer:
            return band.pressure
    return 0.0


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


def collect_cuts(loads, clamp_radius, free_radius):
    cuts = {clamp_radius, free_radius}
    for r_inner, r_outer, _ in loads:
        if clamp_radius < r_inner < free_radius:
            cuts.add(r_inner)
        if clamp_radius < r_outer < free_radius:
            cuts.add(r_outer)
    return sorted(cuts)


def build_regions(loads, clamp_radius, free_radius):
    edges = collect_cuts(loads, clamp_radius, free_radius)
    bands = make_bands(loads)
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


class AnnularPlate:
    def __init__(self, regions, coefficients, rigidity, poisson_ratio, thickness):
        self.regions = regions
        self.coefficients = coefficients
        self.rigidity = rigidity
        self.poisson_ratio = poisson_ratio
        self.thickness = thickness

    def _local(self, r):
        index = region_index(self.regions, r)
        start = index * BASIS_COUNT
        return self.regions[index][2], self.coefficients[start:start + BASIS_COUNT]

    def deflection(self, r):
        pressure, local = self._local(r)
        return float(basis(r) @ local) + particular_deflection(r, pressure, self.rigidity)

    def slope(self, r):
        pressure, local = self._local(r)
        return float(basis_slope(r) @ local) + particular_slope(r, pressure, self.rigidity)

    def curvature(self, r):
        pressure, local = self._local(r)
        return float(basis_curvature(r) @ local) + particular_curvature(r, pressure, self.rigidity)

    def radial_moment(self, r):
        pressure, local = self._local(r)
        homogeneous = float(basis_radial_moment(r, self.rigidity, self.poisson_ratio) @ local)
        return homogeneous + particular_radial_moment(r, pressure, self.poisson_ratio)

    def hoop_moment(self, r):
        pressure, local = self._local(r)
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


def assemble(regions, clamp_radius, free_radius, rigidity, poisson_ratio):
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
        if kind == "curvature":
            return particular_curvature(r, pressure, rigidity)
        if kind == "shear":
            return particular_shear(r, pressure)
        if kind == "radial_moment":
            return particular_radial_moment(r, pressure, poisson_ratio)
        raise ValueError("unknown particular solution")

    matrix[row] = 0.0
    add(0, 1.0, basis(clamp_radius))
    right[row] = -value(0, clamp_radius, "w")
    row += 1

    matrix[row] = 0.0
    add(0, 1.0, basis_slope(clamp_radius))
    right[row] = -value(0, clamp_radius, "slope")
    row += 1

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

    last = len(regions) - 1

    matrix[row] = 0.0
    add(last, 1.0, basis_radial_moment(free_radius, rigidity, poisson_ratio))
    right[row] = -value(last, free_radius, "radial_moment")
    row += 1

    matrix[row] = 0.0
    add(last, 1.0, basis_shear(free_radius, rigidity))
    right[row] = -value(last, free_radius, "shear")

    return matrix, right


def solve(loads, clamp_radius=None, free_radius=None, thickness=None):
    if clamp_radius is None:
        clamp_radius = geometry.BRACKET_RADIUS
    if free_radius is None:
        free_radius = geometry.OUTER_RADIUS
    if thickness is None:
        thickness = geometry.THICKNESS
    rigidity = material.flexural_rigidity(thickness)
    poisson_ratio = material.POISSON_RATIO
    regions = build_regions(loads, clamp_radius, free_radius)
    matrix, right = assemble(regions, clamp_radius, free_radius, rigidity, poisson_ratio)
    coefficients = np.linalg.solve(matrix, right)
    return AnnularPlate(regions, coefficients, rigidity, poisson_ratio, thickness)


def sample(plate, clamp_radius=None, free_radius=None, points=401):
    if clamp_radius is None:
        clamp_radius = geometry.BRACKET_RADIUS
    if free_radius is None:
        free_radius = geometry.OUTER_RADIUS
    radii = np.linspace(clamp_radius, free_radius, points)
    radial_stress = np.array([plate.radial_stress(r) for r in radii])
    hoop_stress = np.array([plate.hoop_stress(r) for r in radii])
    von_mises = np.array([plate.von_mises_stress(r) for r in radii])
    deflection = np.array([plate.deflection(r) for r in radii])
    return radii, radial_stress, hoop_stress, von_mises, deflection


def report(plate, clamp_radius=None, free_radius=None, points=401):
    radii, radial_stress, hoop_stress, von_mises, deflection = sample(
        plate, clamp_radius, free_radius, points
    )
    index = int(np.argmax(von_mises))
    radial_index = int(np.argmax(np.abs(radial_stress)))
    hoop_index = int(np.argmax(np.abs(hoop_stress)))
    deflection_index = int(np.argmax(np.abs(deflection)))
    return {
        "peak_von_mises_stress": float(von_mises[index]),
        "peak_von_mises_radius": float(radii[index]),
        "peak_radial_stress": float(radial_stress[radial_index]),
        "peak_radial_radius": float(radii[radial_index]),
        "peak_hoop_stress": float(hoop_stress[hoop_index]),
        "peak_hoop_radius": float(radii[hoop_index]),
        "max_deflection": float(deflection[deflection_index]),
        "max_deflection_radius": float(radii[deflection_index]),
        "clamp_radial_stress": float(radial_stress[0]),
        "clamp_hoop_stress": float(hoop_stress[0]),
        "free_edge_radial_stress": float(radial_stress[-1]),
        "free_edge_hoop_stress": float(hoop_stress[-1]),
    }


def format_report(summary):
    lines = []
    for key, value in summary.items():
        if key.endswith("_radius"):
            lines.append((key, value, "mm", 1.0e3))
        elif key.endswith("_stress"):
            lines.append((key, value, "MPa", 1.0e-6))
        else:
            lines.append((key, value, "um", 1.0e6))
    width = max(len(name) for name, _, _, _ in lines)
    output = []
    for name, value, unit, scale in lines:
        output.append(name.ljust(width) + "  " + format(value * scale, ".4f").rjust(14) + "  " + unit)
    return "\n".join(output)
