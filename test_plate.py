import math

import numpy as np
import pytest

import beam
import geometry
import loads
import material
import plate


def test_ring_width_matches_diameter_difference():
    assert geometry.RING_WIDTH == pytest.approx((geometry.OUTER_DIAMETER - geometry.INNER_DIAMETER) / 2.0)


def test_overhang_matches_bracket_difference():
    assert geometry.OVERHANG == pytest.approx(geometry.OUTER_RADIUS - geometry.BRACKET_RADIUS)


def test_load_band_inside_plate():
    assert loads.EPOXY_BAND_INNER_RADIUS >= geometry.INNER_RADIUS
    assert loads.EPOXY_BAND_OUTER_RADIUS <= geometry.OUTER_RADIUS


def test_load_band_starts_at_clamp_radius():
    assert loads.EPOXY_BAND_INNER_RADIUS == pytest.approx(geometry.BRACKET_RADIUS)


def test_uniform_pressure_reproduces_total_force():
    pressure = loads.uniform_pressure(loads.PEAK_THRUST)
    total = pressure * loads.bearing_area()
    assert total == pytest.approx(loads.PEAK_THRUST)


def test_ultimate_load_is_safety_factor_times_peak():
    assert loads.ULTIMATE_THRUST == pytest.approx(loads.PEAK_THRUST * loads.SAFETY_FACTOR_ULTIMATE)


def test_zero_load_gives_zero_everything():
    solution = plate.solve([], clamp_radius=geometry.BRACKET_RADIUS, free_radius=geometry.OUTER_RADIUS)
    assert solution.deflection(geometry.OUTER_RADIUS) == pytest.approx(0.0, abs=1e-18)
    assert solution.radial_moment(geometry.BRACKET_RADIUS) == pytest.approx(0.0, abs=1e-9)
    assert solution.hoop_moment(geometry.OUTER_RADIUS) == pytest.approx(0.0, abs=1e-9)


def test_clamped_edge_conditions():
    solution = plate.solve(loads.peak_thrust_load())
    assert solution.deflection(geometry.BRACKET_RADIUS) == pytest.approx(0.0, abs=1e-12)
    assert solution.slope(geometry.BRACKET_RADIUS) == pytest.approx(0.0, abs=1e-9)


def test_free_edge_moment_is_zero():
    solution = plate.solve(loads.peak_thrust_load())
    assert solution.radial_moment(geometry.OUTER_RADIUS) == pytest.approx(0.0, abs=1e-6)


def test_free_edge_shear_is_zero():
    solution = plate.solve(loads.peak_thrust_load())
    step = 1.0e-9
    radius = geometry.OUTER_RADIUS - 1.0e-8

    def shear(r):
        combined = lambda x: solution.curvature(x) + solution.slope(x) / x
        return -solution.rigidity * (combined(r + step) - combined(r - step)) / (2.0 * step)

    assert shear(radius) == pytest.approx(0.0, abs=50.0)


def test_clamp_reaction_equals_applied_load():
    solution = plate.solve(loads.peak_thrust_load())
    step = 1.0e-9
    radius = geometry.BRACKET_RADIUS + 1.0e-8

    def shear(r):
        combined = lambda x: solution.curvature(x) + solution.slope(x) / x
        return -solution.rigidity * (combined(r + step) - combined(r - step)) / (2.0 * step)

    reaction = 2.0 * math.pi * geometry.BRACKET_RADIUS * shear(radius)
    assert reaction == pytest.approx(loads.PEAK_THRUST, rel=1e-5)


def test_two_region_partition_matches_single_region_when_load_is_uniform():
    radius = 0.055
    clamp = 0.035
    pressure = 1.0e5
    whole = plate.solve([(clamp, radius, pressure)], clamp_radius=clamp, free_radius=radius)
    split = plate.solve(
        [(clamp, 0.045, pressure), (0.045, radius, pressure)],
        clamp_radius=clamp,
        free_radius=radius,
    )
    midpoint = 0.050
    assert split.radial_moment(midpoint) == pytest.approx(whole.radial_moment(midpoint), rel=1e-9)
    assert split.hoop_moment(midpoint) == pytest.approx(whole.hoop_moment(midpoint), rel=1e-9)
    assert split.deflection(radius) == pytest.approx(whole.deflection(radius), rel=1e-9)


def test_annulus_uniform_load_matches_symbolic_ground_truth():
    radius = 0.055
    clamp = 0.035
    pressure = 1.0e5
    solution = plate.solve([(clamp, radius, pressure)], clamp_radius=clamp, free_radius=radius)
    assert solution.radial_moment(clamp) == pytest.approx(-25.0060247318874, rel=1e-9)
    assert solution.hoop_moment(clamp) == pytest.approx(-8.251988161522812, rel=1e-9)


def test_peak_stress_is_at_the_clamped_edge():
    solution = plate.solve(loads.peak_thrust_load())
    summary = plate.report(solution)
    assert summary["peak_von_mises_radius"] == pytest.approx(geometry.BRACKET_RADIUS)


def test_peak_stress_is_positive_and_finite():
    solution = plate.solve(loads.peak_thrust_load())
    summary = plate.report(solution)
    assert summary["peak_von_mises_stress"] > 0.0
    assert math.isfinite(summary["peak_von_mises_stress"])


def test_stress_scales_linearly_with_load():
    single = plate.solve(loads.load_bands(loads.PEAK_THRUST))
    double = plate.solve(loads.load_bands(2.0 * loads.PEAK_THRUST))
    assert double.radial_moment(geometry.BRACKET_RADIUS) == pytest.approx(
        2.0 * single.radial_moment(geometry.BRACKET_RADIUS), rel=1e-9
    )


def test_deflection_scales_linearly_with_load():
    single = plate.solve(loads.load_bands(loads.PEAK_THRUST))
    double = plate.solve(loads.load_bands(2.0 * loads.PEAK_THRUST))
    assert double.deflection(geometry.OUTER_RADIUS) == pytest.approx(
        2.0 * single.deflection(geometry.OUTER_RADIUS), rel=1e-9
    )


def test_deflection_scales_inversely_with_cube_of_thickness():
    thin = plate.solve(loads.peak_thrust_load(), thickness=0.006)
    thick = plate.solve(loads.peak_thrust_load(), thickness=0.012)
    ratio = thin.deflection(geometry.OUTER_RADIUS) / thick.deflection(geometry.OUTER_RADIUS)
    assert ratio == pytest.approx(8.0, rel=1e-6)


def test_more_points_does_not_change_peak_stress():
    solution = plate.solve(loads.peak_thrust_load())
    coarse = plate.report(solution, points=51)["peak_von_mises_stress"]
    fine = plate.report(solution, points=2001)["peak_von_mises_stress"]
    assert fine == pytest.approx(coarse, rel=1e-3)


def test_beam_baseline_disagrees_with_plate_by_large_factor():
    solution = plate.solve(loads.peak_thrust_load())
    plate_stress = plate.report(solution)["peak_von_mises_stress"]
    beam_stress = beam.root_bending_stress(loads.PEAK_THRUST)
    assert beam_stress / plate_stress > 10.0


def test_beam_span_is_the_overhang_not_the_ring_width():
    assert beam.span() == pytest.approx(geometry.OVERHANG)
    assert beam.span() != pytest.approx(geometry.RING_WIDTH)


def test_material_flexural_rigidity_positive():
    assert material.flexural_rigidity(geometry.THICKNESS) > 0.0


def test_safety_factor_against_yield_for_peak_load():
    solution = plate.solve(loads.peak_thrust_load())
    peak = plate.report(solution)["peak_von_mises_stress"]
    factor = material.YIELD_STRESS / peak
    assert factor > 1.0
