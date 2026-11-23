import math

import geometry
import loads
import material
import plate


def stress_at_thickness(thickness, force):
    solution = plate.solve(loads.load_bands(force), thickness=thickness)
    summary = plate.report(solution)
    return summary["peak_von_mises_stress"], summary["max_deflection"]


def required_thickness(target_safety_factor, force, reference_thickness=geometry.THICKNESS):
    target_stress = material.YIELD_STRESS / target_safety_factor
    reference_stress, _ = stress_at_thickness(reference_thickness, force)
    return reference_thickness * math.sqrt(reference_stress / target_stress)


def uniform_solution(target_safety_factor, force=loads.ULTIMATE_THRUST, reference_thickness=geometry.THICKNESS):
    thickness = required_thickness(target_safety_factor, force, reference_thickness)
    peak_stress, deflection = stress_at_thickness(thickness, force)
    return {
        "thickness": thickness,
        "peak_stress": peak_stress,
        "safety_factor": material.YIELD_STRESS / peak_stress,
        "deflection": deflection,
        "mass": annulus_mass(thickness),
        "force": force,
    }


def annulus_mass(thickness):
    return material.DENSITY * geometry.plate_volume(thickness=thickness)


def thin(force=loads.ULTIMATE_THRUST, target_safety_factor=loads.SAFETY_FACTOR_ULTIMATE):
    baseline_mass = annulus_mass(geometry.THICKNESS)
    result = uniform_solution(target_safety_factor, force)
    result["baseline_mass"] = baseline_mass
    result["mass_saving"] = baseline_mass - result["mass"]
    result["mass_saving_fraction"] = (baseline_mass - result["mass"]) / baseline_mass
    return result


def report(result):
    lines = []
    lines.append("Optimisation target: safety factor {:.2f} at {:.0f} N".format(
        loads.SAFETY_FACTOR_ULTIMATE, result["force"]
    ))
    lines.append("")
    lines.append("  required thickness      {:10.4f}  mm".format(result["thickness"] * 1e3))
    lines.append("  peak von Mises stress   {:10.4f}  MPa".format(result["peak_stress"] / 1e6))
    lines.append("  achieved safety factor  {:10.4f}".format(result["safety_factor"]))
    lines.append("  maximum deflection      {:10.4f}  um".format(result["deflection"] * 1e6))
    lines.append("")
    lines.append("  baseline mass, 6 mm     {:10.4f}  g".format(result["baseline_mass"] * 1e3))
    lines.append("  optimised mass          {:10.4f}  g".format(result["mass"] * 1e3))
    lines.append("  mass saved              {:10.4f}  g".format(result["mass_saving"] * 1e3))
    lines.append("  mass saving             {:10.4f}  %".format(result["mass_saving_fraction"] * 100.0))
    return "\n".join(lines)
