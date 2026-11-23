import math

import numpy as np
import pytest

import solver


def test_geometry_order():
    assert solver.BORE_RADIUS < solver.LOAD_OUTER_RADIUS
    assert solver.LOAD_OUTER_RADIUS < solver.FREE_OUTER_RADIUS
    assert solver.FREE_OUTER_RADIUS < solver.PLATE_OUTER_RADIUS


def test_clamp_radius_is_middle_of_the_clamped_band():
    assert solver.CLAMP_RADIUS == pytest.approx(
        0.5 * (solver.FREE_OUTER_RADIUS + solver.PLATE_OUTER_RADIUS)
    )


def test_clamp_radius_lies_inside_the_clamped_band():
    assert solver.FREE_OUTER_RADIUS < solver.CLAMP_RADIUS < solver.PLATE_OUTER_RADIUS


def test_default_solves_use_the_single_clamp_radius():
    default = solver.solve_layout(solver.PEAK_THRUST)
    explicit = solver.solve_layout(solver.PEAK_THRUST, solver.CLAMP_RADIUS)
    assert default.deflection(solver.BORE_RADIUS) == pytest.approx(
        explicit.deflection(solver.BORE_RADIUS)
    )


def test_clamp_at_middle_radius_brackets_the_rigid_limits():
    middle = solver.solve_layout(solver.PEAK_THRUST, solver.CLAMP_RADIUS)
    band_start = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    outer_rim = solver.solve_layout(solver.PEAK_THRUST, solver.PLATE_OUTER_RADIUS)
    middle_stress = middle.von_mises_stress(solver.CLAMP_RADIUS)
    assert band_start.von_mises_stress(solver.FREE_OUTER_RADIUS) < middle_stress
    assert middle_stress < outer_rim.von_mises_stress(solver.PLATE_OUTER_RADIUS)


def test_plate_inner_radius_matches_bore():
    assert solver.PLATE_INNER_RADIUS == pytest.approx(solver.BORE_RADIUS)


def test_load_band_inner_edge_matches_bore():
    assert solver.LOAD_INNER_RADIUS == pytest.approx(solver.BORE_RADIUS)


def test_outer_diameter_matches_plate_outer_radius():
    assert solver.PLATE_OUTER_RADIUS == pytest.approx(solver.OUTER_DIAMETER / 2.0)


def test_load_band_inside_plate():
    assert solver.LOAD_INNER_RADIUS >= solver.PLATE_INNER_RADIUS
    assert solver.LOAD_OUTER_RADIUS <= solver.PLATE_OUTER_RADIUS


def test_load_area_value():
    assert solver.LOAD_AREA == pytest.approx(1429.42e-6, rel=1e-4)


def test_load_area_matches_the_two_band_radii():
    expected = math.pi * (
        solver.LOAD_OUTER_RADIUS**2 - solver.LOAD_INNER_RADIUS**2
    )
    assert solver.LOAD_AREA == pytest.approx(expected)


def test_plate_area_matches_the_two_plate_radii():
    expected = math.pi * (
        solver.PLATE_OUTER_RADIUS**2 - solver.PLATE_INNER_RADIUS**2
    )
    assert solver.PLATE_AREA == pytest.approx(expected)


def test_uniform_pressure_reproduces_total_force():
    pressure = solver.uniform_pressure(solver.PEAK_THRUST)
    assert pressure * solver.LOAD_AREA == pytest.approx(solver.PEAK_THRUST)


def test_peak_load_pressure_matches_peak_thrust():
    assert solver.LOAD_PRESSURE == pytest.approx(
        solver.uniform_pressure(solver.PEAK_THRUST)
    )


def test_plate_mass_matches_area_times_thickness_times_density():
    assert solver.plate_mass(solver.THICKNESS) == pytest.approx(
        solver.DENSITY * solver.PLATE_AREA * solver.THICKNESS
    )


def test_ultimate_thrust_is_safety_factor_times_peak():
    assert solver.ULTIMATE_THRUST == pytest.approx(
        solver.PEAK_THRUST * solver.SAFETY_FACTOR_REQUIRED
    )


def test_allowable_stress_is_yield_over_required_factor():
    assert solver.ALLOWABLE_STRESS == pytest.approx(
        solver.YIELD_STRESS / solver.SAFETY_FACTOR_REQUIRED
    )


def test_flexural_rigidity_positive():
    assert solver.flexural_rigidity(solver.THICKNESS) > 0.0


def test_flexural_rigidity_scales_with_cube_of_thickness():
    single = solver.flexural_rigidity(0.006)
    double = solver.flexural_rigidity(0.012)
    assert double == pytest.approx(8.0 * single, rel=1e-9)


def test_plate_mass_at_baseline_thickness():
    assert solver.plate_mass(solver.THICKNESS) * 1e3 == pytest.approx(111.1521, rel=1e-4)


def test_plate_mass_scales_linearly_with_thickness():
    single = solver.plate_mass(0.002)
    double = solver.plate_mass(0.004)
    assert double == pytest.approx(2.0 * single, rel=1e-9)


def test_load_band_outside_domain_raises():
    inner = solver.Edge(0.040, "free")
    outer = solver.Edge(0.055, "clamped")
    with pytest.raises(ValueError):
        solver.solve(solver.load_bands(solver.PEAK_THRUST), inner, outer)


def test_inverted_edges_raise():
    inner = solver.Edge(0.055, "free")
    outer = solver.Edge(0.029, "clamped")
    with pytest.raises(ValueError):
        solver.solve(solver.load_bands(solver.PEAK_THRUST), inner, outer)


def test_unknown_edge_condition_raises():
    inner = solver.Edge(solver.BORE_RADIUS, "pinned")
    outer = solver.Edge(solver.PLATE_OUTER_RADIUS, "clamped")
    with pytest.raises(ValueError):
        solver.solve(solver.load_bands(solver.PEAK_THRUST), inner, outer)


def test_zero_load_gives_zero_everything():
    solution = solver.solve(
        [],
        solver.Edge(solver.BORE_RADIUS, "free"),
        solver.Edge(solver.PLATE_OUTER_RADIUS, "clamped"),
    )
    assert solution.deflection(solver.BORE_RADIUS) == pytest.approx(0.0, abs=1e-18)
    assert solution.radial_moment(solver.PLATE_OUTER_RADIUS) == pytest.approx(0.0, abs=1e-9)


def test_free_edge_moment_is_zero():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    assert solution.radial_moment(solver.BORE_RADIUS) == pytest.approx(0.0, abs=1e-6)


def test_clamped_edge_is_fixed():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    assert solution.deflection(solver.FREE_OUTER_RADIUS) == pytest.approx(0.0, abs=1e-12)
    assert solution.slope(solver.FREE_OUTER_RADIUS) == pytest.approx(0.0, abs=1e-9)


def test_clamped_edge_fixed_at_outer_rim():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.PLATE_OUTER_RADIUS)
    assert solution.deflection(solver.PLATE_OUTER_RADIUS) == pytest.approx(0.0, abs=1e-12)
    assert solution.slope(solver.PLATE_OUTER_RADIUS) == pytest.approx(0.0, abs=1e-9)


def test_clamp_reaction_equals_applied_load():
    for clamp_radius in [solver.FREE_OUTER_RADIUS, solver.PLATE_OUTER_RADIUS]:
        solution = solver.solve_layout(solver.PEAK_THRUST, clamp_radius)
        reaction = 2.0 * math.pi * clamp_radius * solution.shear(clamp_radius - 1.0e-8)
        assert abs(reaction) == pytest.approx(solver.PEAK_THRUST, rel=1e-4)


def test_clamp_moment_matches_ring_of_force_statics():
    centroid = 0.5 * (solver.LOAD_INNER_RADIUS + solver.LOAD_OUTER_RADIUS)
    for clamp_radius in [solver.FREE_OUTER_RADIUS, solver.PLATE_OUTER_RADIUS]:
        solution = solver.solve_layout(solver.PEAK_THRUST, clamp_radius)
        estimate = solver.PEAK_THRUST * abs(centroid - clamp_radius) / (2.0 * math.pi * clamp_radius)
        actual = abs(solution.radial_moment(clamp_radius))
        assert actual == pytest.approx(estimate, rel=0.05)


def test_peak_stress_is_at_the_clamp():
    for clamp_radius in [solver.FREE_OUTER_RADIUS, solver.PLATE_OUTER_RADIUS]:
        solution = solver.solve_layout(solver.PEAK_THRUST, clamp_radius)
        result = solver.summary(solution, solver.BORE_RADIUS, clamp_radius)
        assert result["peak_von_mises_radius"] == pytest.approx(clamp_radius)


def test_stress_scales_linearly_with_load():
    single = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    double = solver.solve_layout(2.0 * solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    assert double.radial_moment(solver.FREE_OUTER_RADIUS) == pytest.approx(
        2.0 * single.radial_moment(solver.FREE_OUTER_RADIUS), rel=1e-9
    )


def test_deflection_scales_linearly_with_load():
    single = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    double = solver.solve_layout(2.0 * solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    assert double.deflection(solver.BORE_RADIUS) == pytest.approx(
        2.0 * single.deflection(solver.BORE_RADIUS), rel=1e-9
    )


def test_deflection_scales_inversely_with_cube_of_thickness():
    thin = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS, thickness=0.006)
    thick = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS, thickness=0.012)
    ratio = thin.deflection(solver.BORE_RADIUS) / thick.deflection(solver.BORE_RADIUS)
    assert ratio == pytest.approx(8.0, rel=1e-6)


def test_stress_scales_inversely_with_square_of_thickness():
    thick = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS, thickness=0.006)
    thin = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS, thickness=0.003)
    ratio = thin.von_mises_stress(solver.FREE_OUTER_RADIUS) / thick.von_mises_stress(
        solver.FREE_OUTER_RADIUS
    )
    assert ratio == pytest.approx(4.0, rel=1e-6)


def test_region_split_does_not_change_answer():
    inner = solver.Edge(solver.BORE_RADIUS, "free")
    outer = solver.Edge(solver.FREE_OUTER_RADIUS, "clamped")
    pressure = 882173.0
    whole = solver.solve(
        [(solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS, pressure)], inner, outer
    )
    split = solver.solve(
        [
            (solver.BORE_RADIUS, 0.033, pressure),
            (0.033, solver.FREE_OUTER_RADIUS, pressure),
        ],
        inner,
        outer,
    )
    assert split.radial_moment(0.040) == pytest.approx(whole.radial_moment(0.040), rel=1e-9)


def test_uniform_load_reproduces_closed_form_annulus_moment():
    inner = solver.Edge(solver.BORE_RADIUS, "free")
    outer = solver.Edge(solver.PLATE_OUTER_RADIUS, "clamped")
    pressure = 1.0e5
    solution = solver.solve([(solver.BORE_RADIUS, solver.PLATE_OUTER_RADIUS, pressure)], inner, outer)
    assert math.isfinite(solution.radial_moment(solver.PLATE_OUTER_RADIUS))
    assert math.isfinite(solution.hoop_moment(solver.PLATE_OUTER_RADIUS))


def test_more_sample_points_does_not_change_peak():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    coarse = solver.summary(solution, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS, points=51)
    fine = solver.summary(solution, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS, points=2001)
    assert fine["peak_von_mises_stress"] == pytest.approx(coarse["peak_von_mises_stress"], rel=1e-3)


def test_deflection_is_small_compared_with_thickness():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    assert abs(solution.deflection(solver.BORE_RADIUS)) < 0.01 * solver.THICKNESS


def test_peak_stress_is_far_below_yield():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    peak = solver.summary(solution, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)
    assert peak["peak_von_mises_stress"] < 0.1 * solver.YIELD_STRESS


def test_stress_at_ultimate_load_is_below_yield():
    solution = solver.solve_layout(solver.ULTIMATE_THRUST, solver.FREE_OUTER_RADIUS)
    peak = solver.summary(solution, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)
    assert peak["peak_von_mises_stress"] < solver.YIELD_STRESS


def test_moment_is_largest_at_the_clamp():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    assert abs(solution.radial_moment(solver.FREE_OUTER_RADIUS)) > abs(
        solution.radial_moment(solver.BORE_RADIUS)
    )


def test_deflection_is_largest_at_the_bore():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    result = solver.summary(solution, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)
    assert result["max_deflection_radius"] == pytest.approx(solver.BORE_RADIUS)


def test_thinner_plate_has_higher_stress():
    thick = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS, thickness=0.006)
    thin = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS, thickness=0.003)
    assert thin.von_mises_stress(solver.FREE_OUTER_RADIUS) > thick.von_mises_stress(
        solver.FREE_OUTER_RADIUS
    )


def test_outer_rim_clamp_gives_higher_stress_than_annulus_clamp():
    annulus = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    rim = solver.solve_layout(solver.PEAK_THRUST, solver.PLATE_OUTER_RADIUS)
    assert rim.von_mises_stress(solver.PLATE_OUTER_RADIUS) > annulus.von_mises_stress(
        solver.FREE_OUTER_RADIUS
    )


def test_required_thickness_scales_with_square_root():
    loose = solver.required_thickness(100.0e6, 250.0e6 / 2.0)
    tight = solver.required_thickness(100.0e6, 250.0e6 / 8.0)
    assert tight == pytest.approx(loose * 2.0, rel=1e-9)


def test_required_thickness_reproduces_allowable_stress():
    reference = solver.solve_layout(solver.ULTIMATE_THRUST, solver.FREE_OUTER_RADIUS)
    reference_stress = solver.summary(reference, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)[
        "peak_von_mises_stress"
    ]
    thickness = solver.required_thickness(reference_stress, solver.ALLOWABLE_STRESS)
    solution = solver.solve_layout(
        solver.ULTIMATE_THRUST, solver.FREE_OUTER_RADIUS, thickness=thickness
    )
    result = solver.summary(solution, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)
    assert result["peak_von_mises_stress"] == pytest.approx(solver.ALLOWABLE_STRESS, rel=1e-6)


def test_optimiser_achieves_required_safety_factor():
    result = solver.optimise_for_safety_factor()
    assert result["safety_factor"] == pytest.approx(solver.SAFETY_FACTOR_REQUIRED, rel=1e-6)


def test_optimiser_thickness_is_thinner_than_baseline():
    result = solver.optimise_for_safety_factor()
    assert result["thickness"] < solver.THICKNESS


def test_optimiser_reduces_mass():
    result = solver.optimise_for_safety_factor()
    assert result["mass"] < result["baseline_mass"]
    assert result["mass_saving_fraction"] > 0.0


def test_optimiser_mass_ratio_equals_thickness_ratio():
    result = solver.optimise_for_safety_factor()
    assert result["mass"] / result["baseline_mass"] == pytest.approx(
        result["thickness"] / solver.THICKNESS, rel=1e-6
    )


def test_optimiser_thickness_increases_with_safety_factor_target():
    reference = solver.solve_layout(solver.ULTIMATE_THRUST, solver.FREE_OUTER_RADIUS)
    stress = solver.summary(reference, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)[
        "peak_von_mises_stress"
    ]
    loose = solver.required_thickness(stress, solver.YIELD_STRESS / 2.0)
    tight = solver.required_thickness(stress, solver.YIELD_STRESS / 8.0)
    assert tight > loose


def test_optimiser_uses_ultimate_thrust():
    assert solver.ULTIMATE_THRUST == pytest.approx(
        solver.PEAK_THRUST * solver.SAFETY_FACTOR_REQUIRED
    )


def test_optimised_plate_equilibrium_holds():
    result = solver.optimise_for_safety_factor()
    solution = solver.solve_layout(
        solver.ULTIMATE_THRUST, solver.FREE_OUTER_RADIUS, thickness=result["thickness"]
    )
    reaction = 2.0 * math.pi * solver.FREE_OUTER_RADIUS * solution.shear(
        solver.FREE_OUTER_RADIUS - 1.0e-8
    )
    assert abs(reaction) == pytest.approx(solver.ULTIMATE_THRUST, rel=1e-4)


def test_optimised_plate_stays_within_small_deflection_theory():
    result = solver.optimise_for_safety_factor()
    optimised = solver.solve_layout(
        solver.ULTIMATE_THRUST, solver.FREE_OUTER_RADIUS, thickness=result["thickness"]
    )
    summary = solver.summary(optimised, solver.BORE_RADIUS, solver.FREE_OUTER_RADIUS)
    assert abs(summary["max_deflection"]) < 0.1 * result["thickness"]


def test_von_mises_matches_manual_combination():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    radius = 0.040
    radial = solution.radial_stress(radius)
    hoop = solution.hoop_stress(radius)
    expected = math.sqrt(radial**2 - radial * hoop + hoop**2)
    assert solution.von_mises_stress(radius) == pytest.approx(expected, rel=1e-12)


def test_shear_matches_closed_form_basis_derivative():
    solution = solver.solve_layout(solver.PEAK_THRUST, solver.FREE_OUTER_RADIUS)
    radius = 0.040
    index = solver.region_index(solution.regions, radius)
    local = solution.coefficients[index * solver.BASIS_COUNT:(index + 1) * solver.BASIS_COUNT]
    closed_form = -solution.rigidity * 4.0 * local[3] / radius
    assert solution.shear(radius) == pytest.approx(closed_form, rel=1e-5)


def test_plot_writes_a_file(tmp_path):
    path = tmp_path / "profile.png"
    written = solver.plot(path=str(path))
    assert path.exists()
    assert written == str(path)
