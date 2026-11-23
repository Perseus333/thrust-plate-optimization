import math

import geometry
import material
import plate_corrected


PEAK_THRUST = 1261.0
SAFETY_FACTOR = 2.5
ULTIMATE_THRUST = PEAK_THRUST * SAFETY_FACTOR

BORE_RADIUS = 0.029
LOAD_OUTER_RADIUS = 0.036
CLAMP_INNER_RADIUS = 0.043
OUTER_RADIUS = 0.055


def load_bands(total_force):
    area = math.pi * (LOAD_OUTER_RADIUS**2 - BORE_RADIUS**2)
    pressure = total_force / area
    return [(BORE_RADIUS, LOAD_OUTER_RADIUS, pressure)]


def solve_layout(total_force, clamp_at_clamp_inner=True, thickness=None):
    bands = load_bands(total_force)
    inner = plate_corrected.Edge(BORE_RADIUS, "free")
    if clamp_at_clamp_inner:
        outer = plate_corrected.Edge(CLAMP_INNER_RADIUS, "clamped")
    else:
        outer = plate_corrected.Edge(OUTER_RADIUS, "clamped")
    return plate_corrected.solve(bands, inner, outer, thickness=thickness)


def check_edge_conditions(solution, r_inner, r_outer, inner_condition, outer_condition):
    lines = []
    if inner_condition == "free":
        lines.append(("Mr at inner edge", solution.radial_moment(r_inner)))
    if outer_condition == "clamped":
        lines.append(("w at outer edge", solution.deflection(r_outer)))
        lines.append(("slope at outer edge", solution.slope(r_outer)))
    return lines
