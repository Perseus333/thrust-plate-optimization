import math

import numpy as np
import pytest

import solver


def field_peak(force, kind, clamp=solver.CLAMP_RADIUS):
    solution = solver.solve(total_force=force, clamp_radius=clamp)
    return solver.peak(solution, kind, solver.BORE_RADIUS, clamp)


def test_geometry_order():
    assert solver.BORE_RADIUS < solver.LOAD_OUTER_RADIUS
    assert solver.LOAD_OUTER_RADIUS < solver.CLAMP_RADIUS
    assert solver.CLAMP_RADIUS < solver.PLATE_OUTER_RADIUS


def test_clamp_radius_is_middle_of_the_clamped_band():
    assert solver.CLAMP_RADIUS == pytest.approx(0.5 * (0.043 + 0.055))


def test_clamp_radius_lies_inside_the_clamped_band():
    assert 0.043 < solver.CLAMP_RADIUS < solver.PLATE_OUTER_RADIUS


def test_outer_diameter_matches_plate_outer_radius():
    assert solver.PLATE_OUTER_RADIUS == pytest.approx(solver.OUTER_DIAMETER / 2.0)


def test_load_area_value():
    assert solver.LOAD_AREA == pytest.approx(1429.42e-6, rel=1e-4)


def test_load_area_matches_the_two_band_radii():
    expected = math.pi * (solver.LOAD_OUTER_RADIUS**2 - solver.BORE_RADIUS**2)
    assert solver.LOAD_AREA == pytest.approx(expected)


def test_plate_area_matches_the_two_plate_radii():
    expected = math.pi * (solver.PLATE_OUTER_RADIUS**2 - solver.BORE_RADIUS**2)
    assert solver.PLATE_AREA == pytest.approx(expected)


def test_peak_load_pressure_matches_peak_thrust():
    assert solver.PEAK_THRUST / solver.LOAD_AREA == pytest.approx(882173.1, rel=1e-5)


def test_ultimate_thrust_is_safety_factor_times_peak():
    assert solver.ULTIMATE_THRUST == pytest.approx(
        solver.PEAK_THRUST * solver.SAFETY_FACTOR_REQUIRED
    )


def test_allowable_stress_is_yield_over_required_factor():
    assert solver.ALLOWABLE_STRESS == pytest.approx(
        solver.YIELD_STRESS / solver.SAFETY_FACTOR_REQUIRED
    )


def test_stiffness_scales_with_cube_of_thickness():
    assert solver.stiffness(0.012) == pytest.approx(8.0 * solver.stiffness(0.006), rel=1e-9)


def test_default_solves_use_the_single_clamp_radius():
    default = solver.solve()
    explicit = solver.solve(clamp_radius=solver.CLAMP_RADIUS)
    assert solver.field(default, solver.BORE_RADIUS, "deflection") == pytest.approx(
        solver.field(explicit, solver.BORE_RADIUS, "deflection")
    )


def test_clamp_at_middle_radius_brackets_the_rigid_limits():
    middle = field_peak(solver.PEAK_THRUST, "von_mises", solver.CLAMP_RADIUS)[1]
    band_start = field_peak(solver.PEAK_THRUST, "von_mises", 0.043)[1]
    outer_rim = field_peak(solver.PEAK_THRUST, "von_mises", 0.055)[1]
    assert band_start < middle < outer_rim


def test_free_edge_moment_is_zero():
    solution = solver.solve(clamp_radius=0.043)
    assert solver.field(solution, solver.BORE_RADIUS, "radial_moment") == pytest.approx(
        0.0, abs=1e-6
    )


def test_clamped_edge_is_fixed():
    for clamp in (0.043, 0.055):
        solution = solver.solve(clamp_radius=clamp)
        assert solver.field(solution, clamp, "deflection") == pytest.approx(0.0, abs=1e-12)
        assert solver.field(solution, clamp, "slope") == pytest.approx(0.0, abs=1e-9)


def test_clamp_reaction_equals_applied_load():
    for clamp in (0.043, 0.055):
        solution = solver.solve(clamp_radius=clamp)
        reaction = 2.0 * math.pi * clamp * solver.field(
            solution, clamp - 1.0e-8, "shear"
        )
        assert abs(reaction) == pytest.approx(solver.PEAK_THRUST, rel=1e-4)


def test_clamp_moment_matches_ring_of_force_statics():
    centroid = 0.5 * (solver.BORE_RADIUS + solver.LOAD_OUTER_RADIUS)
    for clamp in (0.043, 0.055):
        solution = solver.solve(clamp_radius=clamp)
        estimate = solver.PEAK_THRUST * abs(centroid - clamp) / (2.0 * math.pi * clamp)
        actual = abs(solver.field(solution, clamp, "radial_moment"))
        assert actual == pytest.approx(estimate, rel=0.05)


def test_peak_stress_is_at_the_clamp():
    for clamp in (0.043, 0.055):
        radius = field_peak(solver.PEAK_THRUST, "von_mises", clamp)[0]
        assert radius == pytest.approx(clamp)


def test_stress_scales_linearly_with_load():
    single = solver.solve(clamp_radius=0.043)
    double = solver.solve(clamp_radius=0.043, total_force=2.0 * solver.PEAK_THRUST)
    assert solver.field(double, 0.043, "radial_moment") == pytest.approx(
        2.0 * solver.field(single, 0.043, "radial_moment"), rel=1e-9
    )


def test_deflection_scales_linearly_with_load():
    single = solver.solve(clamp_radius=0.043)
    double = solver.solve(clamp_radius=0.043, total_force=2.0 * solver.PEAK_THRUST)
    assert solver.field(double, solver.BORE_RADIUS, "deflection") == pytest.approx(
        2.0 * solver.field(single, solver.BORE_RADIUS, "deflection"), rel=1e-9
    )


def test_deflection_scales_inversely_with_cube_of_thickness():
    thin = solver.solve(clamp_radius=0.043, thickness=0.006)
    thick = solver.solve(clamp_radius=0.043, thickness=0.012)
    ratio = solver.field(thin, solver.BORE_RADIUS, "deflection") / solver.field(
        thick, solver.BORE_RADIUS, "deflection"
    )
    assert ratio == pytest.approx(8.0, rel=1e-6)


def test_stress_scales_inversely_with_square_of_thickness():
    thick = solver.solve(clamp_radius=0.043, thickness=0.006)
    thin = solver.solve(clamp_radius=0.043, thickness=0.003)
    ratio = solver.field(thin, 0.043, "von_mises") / solver.field(thick, 0.043, "von_mises")
    assert ratio == pytest.approx(4.0, rel=1e-6)


def test_more_sample_points_does_not_change_peak():
    solution = solver.solve(clamp_radius=0.043)
    coarse = solver.peak(solution, "von_mises", solver.BORE_RADIUS, 0.043)[1]
    radii = np.linspace(solver.BORE_RADIUS, 0.043, 2001)
    fine = max(solver.field(solution, r, "von_mises") for r in radii)
    assert fine == pytest.approx(coarse, rel=1e-3)


def test_deflection_is_small_compared_with_thickness():
    solution = solver.solve(clamp_radius=0.043)
    assert abs(solver.field(solution, solver.BORE_RADIUS, "deflection")) < (
        0.01 * solver.THICKNESS
    )


def test_peak_stress_is_far_below_yield():
    assert field_peak(solver.PEAK_THRUST, "von_mises", 0.043)[1] < 0.1 * solver.YIELD_STRESS


def test_stress_at_ultimate_load_is_below_yield():
    assert field_peak(solver.ULTIMATE_THRUST, "von_mises", 0.043)[1] < solver.YIELD_STRESS


def test_moment_is_largest_at_the_clamp():
    solution = solver.solve(clamp_radius=0.043)
    assert abs(solver.field(solution, 0.043, "radial_moment")) > abs(
        solver.field(solution, solver.BORE_RADIUS, "radial_moment")
    )


def test_deflection_is_largest_at_the_bore():
    radius = field_peak(solver.PEAK_THRUST, "deflection", 0.043)[0]
    assert radius == pytest.approx(solver.BORE_RADIUS)


def test_thinner_plate_has_higher_stress():
    thick = solver.solve(clamp_radius=0.043, thickness=0.006)
    thin = solver.solve(clamp_radius=0.043, thickness=0.003)
    assert solver.field(thin, 0.043, "von_mises") > solver.field(thick, 0.043, "von_mises")


def test_outer_rim_clamp_gives_higher_stress_than_annulus_clamp():
    annulus = field_peak(solver.PEAK_THRUST, "von_mises", 0.043)[1]
    rim = field_peak(solver.PEAK_THRUST, "von_mises", 0.055)[1]
    assert rim > annulus


def test_von_mises_matches_manual_combination():
    solution = solver.solve(clamp_radius=0.043)
    radius = 0.040
    radial = solver.field(solution, radius, "radial_stress")
    hoop = solver.field(solution, radius, "hoop_stress")
    expected = math.sqrt(radial**2 - radial * hoop + hoop**2)
    assert solver.field(solution, radius, "von_mises") == pytest.approx(expected, rel=1e-12)


def test_shear_matches_closed_form_basis_derivative():
    solution = solver.solve(clamp_radius=0.043)
    for radius in (0.032, 0.040):
        region = 0 if radius < solver.LOAD_OUTER_RADIUS else 1
        c4 = solution["coefficients"][region * 4 + 3]
        pressure = solution["pressure"] if region == 0 else 0.0
        expected = -solution["rigidity"] * 4.0 * c4 / radius + pressure * (-radius / 2.0)
        assert solver.field(solution, radius, "shear") == pytest.approx(expected, rel=1e-12)


def test_shear_is_zero_at_the_free_bore():
    for clamp in (0.043, 0.055):
        solution = solver.solve(clamp_radius=clamp)
        assert solver.field(solution, solver.BORE_RADIUS, "shear") == pytest.approx(
            0.0, abs=1e-9
        )


def test_zero_load_gives_zero_everything():
    solution = solver.solve(clamp_radius=0.055, total_force=0.0)
    assert solver.field(solution, solver.BORE_RADIUS, "deflection") == pytest.approx(
        0.0, abs=1e-18
    )
    assert solver.field(solution, solver.PLATE_OUTER_RADIUS, "radial_moment") == pytest.approx(
        0.0, abs=1e-9
    )


def test_clamp_outside_the_plate_raises():
    with pytest.raises(ValueError):
        solver.solve(clamp_radius=0.060)


def test_clamp_at_or_inside_the_bore_raises():
    with pytest.raises(ValueError):
        solver.solve(clamp_radius=solver.BORE_RADIUS)


def test_region_boundary_belongs_to_the_outer_region():
    solution = solver.solve(clamp_radius=0.043)
    radius = solver.LOAD_OUTER_RADIUS
    outer = solver.term_at(solution, radius, "deflection")
    inner = solver.term_at(solution, radius - 1.0e-12, "deflection")
    assert outer == pytest.approx(inner, rel=1e-6)


def test_optimise_achieves_required_safety_factor():
    result = solver.optimise()
    assert result["safety_factor"] == pytest.approx(solver.SAFETY_FACTOR_REQUIRED, rel=1e-6)


def test_optimise_thickness_is_thinner_than_baseline():
    assert solver.optimise()["thickness"] < solver.THICKNESS


def test_optimise_reduces_mass():
    result = solver.optimise()
    assert result["mass"] < result["baseline_mass"]
    assert result["mass_saving_fraction"] > 0.0


def test_optimise_mass_ratio_equals_thickness_ratio():
    result = solver.optimise()
    assert result["mass"] / result["baseline_mass"] == pytest.approx(
        result["thickness"] / solver.THICKNESS, rel=1e-6
    )


def test_optimise_uses_ultimate_thrust():
    assert solver.ULTIMATE_THRUST == pytest.approx(
        solver.PEAK_THRUST * solver.SAFETY_FACTOR_REQUIRED
    )


def test_optimised_plate_equilibrium_holds():
    result = solver.optimise()
    solution = solver.solve(
        thickness=result["thickness"], clamp_radius=0.043, total_force=solver.ULTIMATE_THRUST
    )
    reaction = 2.0 * math.pi * 0.043 * solver.field(solution, 0.043 - 1.0e-8, "shear")
    assert abs(reaction) == pytest.approx(solver.ULTIMATE_THRUST, rel=1e-4)


def test_optimised_plate_stays_within_small_deflection_theory():
    result = solver.optimise()
    solution = solver.solve(
        thickness=result["thickness"], clamp_radius=0.043, total_force=solver.ULTIMATE_THRUST
    )
    deflection = solver.peak(solution, "deflection", solver.BORE_RADIUS, 0.043)[1]
    assert abs(deflection) < 0.1 * result["thickness"]


def test_plate_mass_relation_at_baseline():
    assert solver.DENSITY * solver.PLATE_AREA * solver.THICKNESS * 1e3 == pytest.approx(
        111.1521, rel=1e-4
    )


def test_report_prints_the_documented_numbers(capsys):
    solver.report()
    out = capsys.readouterr().out
    assert "9.7796 MPa" in out
    assert "25.56" in out
    assert "2.9668 mm" in out
    assert "50.5541 %" in out


def test_plot_writes_a_file(tmp_path):
    path = tmp_path / "profile.png"
    written = solver.plot(path=str(path))
    assert path.exists()
    assert written == str(path)
