import math

import geometry
import loads
import material


CLAMP_FACE_RADIUS = geometry.BRACKET_RADIUS
FREE_END_RADIUS = geometry.OUTER_RADIUS


def span():
    return FREE_END_RADIUS - CLAMP_FACE_RADIUS


def section_modulus(width=None, thickness=None):
    if width is None:
        width = geometry.RING_WIDTH
    if thickness is None:
        thickness = geometry.THICKNESS
    return width * thickness**2 / 6.0


def tip_deflection(force=None, width=None, thickness=None):
    if force is None:
        force = loads.PEAK_THRUST
    if width is None:
        width = geometry.RING_WIDTH
    if thickness is None:
        thickness = geometry.THICKNESS
    inertia = geometry.second_moment_of_area(width, thickness)
    rigidity = material.YOUNGS_MODULUS * inertia
    return force * span() ** 3 / (3.0 * rigidity)


def root_bending_stress(force=None, width=None, thickness=None):
    if force is None:
        force = loads.PEAK_THRUST
    if width is None:
        width = geometry.RING_WIDTH
    if thickness is None:
        thickness = geometry.THICKNESS
    moment = force * span()
    return 6.0 * moment / (width * thickness**2)


def report(force=None):
    if force is None:
        force = loads.PEAK_THRUST
    stress = root_bending_stress(force)
    return {
        "beam_span": span(),
        "cantilever_tip_deflection": tip_deflection(force),
        "cantilever_root_bending_stress": stress,
        "cantilever_safety_factor": material.YIELD_STRESS / abs(stress),
    }
