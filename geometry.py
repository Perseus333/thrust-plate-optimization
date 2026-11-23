import math

OUTER_DIAMETER = 0.110
INNER_DIAMETER = 0.058
BRACKET_OUTER_DIAMETER = 0.070
THICKNESS = 0.006
CABLE_HOLE_DIAMETER = 0.008

OUTER_RADIUS = OUTER_DIAMETER / 2.0
INNER_RADIUS = INNER_DIAMETER / 2.0
BRACKET_RADIUS = BRACKET_OUTER_DIAMETER / 2.0

RING_WIDTH = OUTER_RADIUS - INNER_RADIUS
OVERHANG = OUTER_RADIUS - BRACKET_RADIUS


def annulus_area(r_inner, r_outer):
    return math.pi * (r_outer**2 - r_inner**2)


def plate_volume(r_outer=OUTER_RADIUS, r_inner=INNER_RADIUS, thickness=THICKNESS):
    return annulus_area(r_inner, r_outer) * thickness


def second_moment_of_area(width, thickness):
    return width * thickness**3 / 12.0
