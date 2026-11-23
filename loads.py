import math

PEAK_THRUST = 1261.0
SAFETY_FACTOR_ULTIMATE = 2.5
ULTIMATE_THRUST = PEAK_THRUST * SAFETY_FACTOR_ULTIMATE

EPOXY_BAND_INNER_RADIUS = 0.035
EPOXY_BAND_OUTER_RADIUS = 0.050

BURN_TIME = 2.7
TOTAL_IMPULSE = 2714.0


def bearing_area(r_inner=EPOXY_BAND_INNER_RADIUS, r_outer=EPOXY_BAND_OUTER_RADIUS):
    return math.pi * (r_outer**2 - r_inner**2)


def uniform_pressure(total_force, r_inner=EPOXY_BAND_INNER_RADIUS, r_outer=EPOXY_BAND_OUTER_RADIUS):
    return total_force / bearing_area(r_inner, r_outer)


def load_bands(total_force, r_inner=EPOXY_BAND_INNER_RADIUS, r_outer=EPOXY_BAND_OUTER_RADIUS):
    pressure = uniform_pressure(total_force, r_inner, r_outer)
    return [(r_inner, r_outer, pressure)]


def peak_thrust_load():
    return load_bands(PEAK_THRUST)


def ultimate_thrust_load():
    return load_bands(ULTIMATE_THRUST)
