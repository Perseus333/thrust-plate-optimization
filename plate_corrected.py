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


def assemble(regions, r_inner, r_outer, inner_condition, outer_condition,
             rigidity, poisson_ratio):
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

    first_edge = regions[0][0]
    last_edge = regions[-1][1]
    add_edge_condition(0, first_edge, inner_condition)
    add_edge_condition(len(regions) - 1, last_edge, outer_condition)

    return matrix, right


class Edge:
    def __init__(self, radius, condition):
        self.radius = radius
        self.condition = condition


def solve(bands, inner_edge, outer_edge, thickness=None):
    if thickness is None:
        thickness = geometry.THICKNESS
    if inner_edge.radius >= outer_edge.radius:
        raise ValueError("inner edge radius must be smaller than outer edge radius")
    for band_inner, band_outer, _ in bands:
        if band_inner < inner_edge.radius or band_outer > outer_edge.radius:
            raise ValueError("load band lies outside the solved domain")
    rigidity = material.flexural_rigidity(thickness)
    poisson_ratio = material.POISSON_RATIO
    regions = build_regions(bands, inner_edge.radius, outer_edge.radius)
    matrix, right = assemble(
        regions,
        inner_edge.radius,
        outer_edge.radius,
        inner_edge.condition,
        outer_edge.condition,
        rigidity,
        poisson_ratio,
    )
    coefficients = np.linalg.solve(matrix, right)
    return AnnularPlate(regions, coefficients, rigidity, poisson_ratio, thickness)


def sample(plate, r_start, r_end, points=401):
    radii = np.linspace(r_start, r_end, points)
    radial_stress = np.array([plate.radial_stress(r) for r in radii])
    hoop_stress = np.array([plate.hoop_stress(r) for r in radii])
    von_mises = np.array([plate.von_mises_stress(r) for r in radii])
    deflection = np.array([plate.deflection(r) for r in radii])
    return radii, radial_stress, hoop_stress, von_mises, deflection


def report(plate, r_start, r_end, points=401):
    radii, radial_stress, hoop_stress, von_mises, deflection = sample(
        plate, r_start, r_end, points
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
    }
