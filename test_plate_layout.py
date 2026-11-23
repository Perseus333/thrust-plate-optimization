import math

import pytest

import geometry
import material
import plate_corrected
import plate_layout


def test_corrected_geometry_order():
    assert plate_layout.BORE_RADIUS < plate_layout.LOAD_OUTER_RADIUS
    assert plate_layout.LOAD_OUTER_RADIUS < plate_layout.CLAMP_INNER_RADIUS
    assert plate_layout.CLAMP_INNER_RADIUS < plate_layout.OUTER_RADIUS


def test_load_annulus_area():
    area = math.pi * (plate_layout.LOAD_OUTER_RADIUS**2 - plate_layout.BORE_RADIUS**2)
    assert area == pytest.approx(1429.42e-6, rel=1e-4)


def test_pressure_reproduces_total_force():
    bands = plate_layout.load_bands(1261.0)
    r_inner, r_outer, pressure = bands[0]
    area = math.pi * (r_outer**2 - r_inner**2)
    assert pressure * area == pytest.approx(1261.0)


def test_load_band_outside_domain_raises():
    inner = plate_corrected.Edge(0.040, "free")
    outer = plate_corrected.Edge(0.055, "clamped")
    bands = plate_layout.load_bands(1261.0)
    with pytest.raises(ValueError):
        plate_corrected.solve(bands, inner, outer)


def test_free_edge_moment_is_zero():
    solution = plate_layout.solve_layout(1261.0)
    assert solution.radial_moment(plate_layout.BORE_RADIUS) == pytest.approx(0.0, abs=1e-6)


def test_clamped_edge_is_fixed():
    solution = plate_layout.solve_layout(1261.0)
    assert solution.deflection(plate_layout.CLAMP_INNER_RADIUS) == pytest.approx(0.0, abs=1e-12)
    assert solution.slope(plate_layout.CLAMP_INNER_RADIUS) == pytest.approx(0.0, abs=1e-9)


def test_clamp_reaction_equals_applied_load():
    solution = plate_layout.solve_layout(1261.0)
    step = 1.0e-9
    radius = plate_layout.CLAMP_INNER_RADIUS - 1.0e-8

    def shear(r):
        combined = lambda x: solution.curvature(x) + solution.slope(x) / x
        return -solution.rigidity * (combined(r + step) - combined(r - step)) / (2.0 * step)

    reaction = 2.0 * math.pi * plate_layout.CLAMP_INNER_RADIUS * shear(radius)
    assert abs(reaction) == pytest.approx(1261.0, rel=1e-5)


def test_clamp_at_outer_rim_equilibrium():
    solution = plate_layout.solve_layout(1261.0, clamp_at_clamp_inner=False)
    step = 1.0e-9
    radius = plate_layout.OUTER_RADIUS - 1.0e-8

    def shear(r):
        combined = lambda x: solution.curvature(x) + solution.slope(x) / x
        return -solution.rigidity * (combined(r + step) - combined(r - step)) / (2.0 * step)

    reaction = 2.0 * math.pi * plate_layout.OUTER_RADIUS * shear(radius)
    assert abs(reaction) == pytest.approx(1261.0, rel=1e-5)


def test_stress_scales_linearly_with_load():
    single = plate_layout.solve_layout(1261.0)
    double = plate_layout.solve_layout(2.0 * 1261.0)
    assert double.radial_moment(plate_layout.CLAMP_INNER_RADIUS) == pytest.approx(
        2.0 * single.radial_moment(plate_layout.CLAMP_INNER_RADIUS), rel=1e-9
    )


def test_deflection_scales_inversely_with_cube_of_thickness():
    thin = plate_layout.solve_layout(1261.0, thickness=0.006)
    thick = plate_layout.solve_layout(1261.0, thickness=0.012)
    ratio = thin.deflection(plate_layout.BORE_RADIUS) / thick.deflection(plate_layout.BORE_RADIUS)
    assert ratio == pytest.approx(8.0, rel=1e-6)


def test_peak_stress_is_at_the_clamp():
    solution = plate_layout.solve_layout(1261.0)
    summary = plate_corrected.report(solution, plate_layout.BORE_RADIUS, plate_layout.CLAMP_INNER_RADIUS)
    assert summary["peak_von_mises_radius"] == pytest.approx(plate_layout.CLAMP_INNER_RADIUS)


def test_region_split_does_not_change_answer():
    inner = plate_corrected.Edge(plate_layout.BORE_RADIUS, "free")
    outer = plate_corrected.Edge(plate_layout.CLAMP_INNER_RADIUS, "clamped")
    pressure = 882173.0
    whole = plate_corrected.solve(
        [(plate_layout.BORE_RADIUS, plate_layout.CLAMP_INNER_RADIUS, pressure)], inner, outer
    )
    split = plate_corrected.solve(
        [(plate_layout.BORE_RADIUS, 0.033, pressure), (0.033, plate_layout.CLAMP_INNER_RADIUS, pressure)],
        inner,
        outer,
    )
    assert split.radial_moment(0.040) == pytest.approx(whole.radial_moment(0.040), rel=1e-9)


def test_thinner_plate_has_higher_stress():
    thick = plate_layout.solve_layout(1261.0, thickness=0.006)
    thin = plate_layout.solve_layout(1261.0, thickness=0.003)
    assert thin.von_mises_stress(plate_layout.CLAMP_INNER_RADIUS) > thick.von_mises_stress(
        plate_layout.CLAMP_INNER_RADIUS
    )


def test_optimised_thickness_reaches_target():
    solution = plate_layout.solve_layout(plate_layout.ULTIMATE_THRUST)
    reference = plate_corrected.report(solution, plate_layout.BORE_RADIUS, plate_layout.CLAMP_INNER_RADIUS)
    target = material.YIELD_STRESS / 2.5
    thickness = 6.0e-3 * math.sqrt(reference["peak_von_mises_stress"] / target)
    optimised = plate_layout.solve_layout(plate_layout.ULTIMATE_THRUST, thickness=thickness)
    summary = plate_corrected.report(optimised, plate_layout.BORE_RADIUS, plate_layout.CLAMP_INNER_RADIUS)
    factor = material.YIELD_STRESS / summary["peak_von_mises_stress"]
    assert factor == pytest.approx(2.5, rel=1e-4)


def test_outer_rim_clamp_gives_higher_stress_than_annulus_clamp():
    annulus = plate_layout.solve_layout(1261.0, clamp_at_clamp_inner=True)
    rim = plate_layout.solve_layout(1261.0, clamp_at_clamp_inner=False)
    assert rim.von_mises_stress(0.050) > annulus.von_mises_stress(0.040)
