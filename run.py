import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import beam
import geometry
import loads
import material
import plate


def build_summary():
    solution = plate.solve(loads.peak_thrust_load())
    summary = plate.report(solution)
    radii, radial_stress, hoop_stress, von_mises, deflection = plate.sample(solution)
    return solution, summary, radii, radial_stress, hoop_stress, von_mises, deflection


def print_summary(summary):
    print("Annular thrust plate, Al 6082-T6, t = 6 mm")
    print("Clamped at r = {:.1f} mm, free at r = {:.1f} mm".format(
        geometry.BRACKET_RADIUS * 1e3, geometry.OUTER_RADIUS * 1e3
    ))
    print("Load {:.0f} N over the epoxy band {:.0f}-{:.0f} mm".format(
        loads.PEAK_THRUST,
        loads.EPOXY_BAND_INNER_RADIUS * 1e3,
        loads.EPOXY_BAND_OUTER_RADIUS * 1e3,
    ))
    print("Equivalent pressure {:.1f} kPa over {:.1f} mm^2".format(
        loads.uniform_pressure(loads.PEAK_THRUST) / 1e3, loads.bearing_area() * 1e6
    ))
    print()
    print(plate.format_report(summary))
    print()
    peak = summary["peak_von_mises_stress"]
    print("Safety factor against yield {:.0f} MPa : {:.2f}".format(
        material.YIELD_STRESS / 1e6, material.YIELD_STRESS / peak
    ))
    print("Safety factor against ultimate {:.0f} MPa : {:.2f}".format(
        material.ULTIMATE_STRESS / 1e6, material.ULTIMATE_STRESS / peak
    ))
    print()
    beam_stress = beam.root_bending_stress(loads.PEAK_THRUST)
    beam_deflection = beam.tip_deflection(loads.PEAK_THRUST)
    print("Clamped strip baseline (span {:.1f} mm):".format(beam.span() * 1e3))
    print("  root bending stress  {:10.2f}  MPa".format(beam_stress / 1e6))
    print("  tip deflection       {:10.2f}  um".format(beam_deflection * 1e6))
    print("  safety factor        {:10.2f}".format(material.YIELD_STRESS / beam_stress))
    print()
    print("Beam over-predicts stress by a factor of {:.1f}".format(
        beam_stress / peak
    ))


def plot(radii, radial_stress, hoop_stress, von_mises, deflection, path):
    millimetres = radii * 1e3
    figure, axes = plt.subplots(2, 1, figsize=(7.0, 8.0))
    axes[0].plot(millimetres, radial_stress / 1e6, label="radial stress")
    axes[0].plot(millimetres, hoop_stress / 1e6, label="hoop stress")
    axes[0].plot(millimetres, von_mises / 1e6, label="von Mises", linestyle="--")
    axes[0].axhline(material.YIELD_STRESS / 1e6, color="red", linestyle=":", label="yield 6082-T6")
    axes[0].axvline(loads.EPOXY_BAND_OUTER_RADIUS * 1e3, color="grey", linestyle="-.", label="load band edge")
    axes[0].set_xlabel("radius, mm")
    axes[0].set_ylabel("bending stress, MPa")
    axes[0].set_title("Through-width bending stress, 1261 N")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].plot(millimetres, deflection * 1e6, color="black")
    axes[1].set_xlabel("radius, mm")
    axes[1].set_ylabel("deflection, um")
    axes[1].set_title("Transverse deflection, 1261 N")
    axes[1].grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main():
    solution, summary, radii, radial_stress, hoop_stress, von_mises, deflection = build_summary()
    print_summary(summary)
    plot(radii, radial_stress, hoop_stress, von_mises, deflection, "stress_profile.png")
    print()
    print("wrote stress_profile.png")


if __name__ == "__main__":
    main()
