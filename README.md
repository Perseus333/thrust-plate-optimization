# Rocket Thrust-plate Analytical Optimisation

Analytical stress and deflection solver for the annular thrust plate of a high-power
rocket, used as the basis for later geometry optimisation.

## The part

A flat annular plate, clamped around its inner region by a bracket and loaded
transversely by the motor thrust carried through an epoxy band.

| Property | Value |
| :---- | :---- |
| Outer diameter | 110 mm |
| Inner diameter | 58 mm |
| Bracket outer diameter | 70 mm |
| Thickness | 6 mm |
| Material | Al 6082-T6 |
| Peak thrust | 1261 N |

Derived geometry:

| Quantity | Value |
| :---- | :---- |
| Radial width of the ring | 26.0 mm |
| Clamp radius | 35.0 mm |
| Free (outer) edge radius | 55.0 mm |
| Unsupported overhang | 20.0 mm |
| Load band (epoxy) | 35.0 - 50.0 mm |
| Bearing area | 4005.5 mm^2 |
| Equivalent pressure at 1261 N | 314.8 kPa |

## Why this is a plate and not a beam

The part is a closed annular ring, not a cantilever strip, so three things differ
from a beam idealisation:

1. The load path runs around the circumference in hoop bending, which a radial
   strip cannot represent.
2. The clamp acts at r = 35 mm, not at the inner bore r = 29 mm.
3. The peak bending moment occurs at the clamped edge, not at the hole.

The clamped-strip baseline in `beam.py` is retained for comparison. It
over-predicts the peak stress by a factor of about 26.

## Method

Axisymmetric Kirchhoff (thin) plate bending. The transverse deflection satisfies

    D * laplacian(laplacian(w)) = q(r),      D = E t^3 / (12 (1 - nu^2))

and the homogeneous solution is

    w(r) = C1 + C2 r^2 + C3 log(r) + C4 r^2 log(r)

with the particular integral for a uniform band pressure

    w_p(r) = q r^4 / (64 D)

The bending moments follow from the curvature relations

    M_r = -D ( w'' + nu w' / r )
    M_t = -D ( w' / r + nu w'' )

and the transverse shear from

    Q = -D d/dr ( w'' + w' / r )

Bending stresses are recovered from

    sigma_r = 6 M_r / t^2
    sigma_t = 6 M_t / t^2

The domain is split into regions at every load-band discontinuity. Within each
region the solution is the homogeneous form plus its own particular integral, and
the regions are joined by continuity of w, w', M_r and Q. The four remaining
conditions are

- clamped edge: w = 0 and w' = 0 at r = 35 mm
- free edge: M_r = 0 and Q = 0 at r = 55 mm

This gives a small dense linear system, solved directly. The radial derivative
needed for the shear condition is available in closed form, since for the
homogeneous basis

    d/dr ( w'' + w' / r ) = 4 C4 / r

so no numerical differentiation is used anywhere in the solver.

Note that Sadd, *Elasticity*, 3rd edition, contains no plate bending theory
(its two-dimensional solutions are plane stress and plane strain, which are
in-plane). The plate formulation above is therefore not sourced from that text.
The relevant overlap with the book is the Lame-type two-constant boundary
condition structure of Section 8.4.1 / Example 8.6, and the polar coordinate
equilibrium and strain relations of Sections 7.6 and 3.8.

## Results at 1261 N

| Quantity | Value |
| :---- | :---- |
| Peak von Mises stress | 6.19 MPa at r = 35 mm |
| Peak radial stress | -7.01 MPa |
| Peak hoop stress | -2.31 MPa |
| Maximum deflection | 2.10 um at r = 55 mm |
| Safety factor against yield (250 MPa) | 40.4 |
| Safety factor against ultimate (300 MPa) | 48.5 |

At the 2.5x ultimate load of 3152 N the peak stress is 15.47 MPa, a safety factor
of 16.2 against yield.

The peak deflection is 0.035 percent of the plate thickness, so the small
deflection and linearity assumptions are satisfied with a very large margin.

## Validation

No ANSYS results are available yet. The following independent checks are in
`test_plate.py`, all passing:

| Check | Agreement |
| :---- | :---- |
| Exact symbolic annulus solution, uniform load (Mr, Mt at clamp) | 1e-9 relative |
| Global equilibrium: clamp shear reaction equals applied load | 1e-5 relative |
| Free edge conditions Mr = 0 and Q = 0 | exact to roundoff |
| Clamped edge conditions w = 0 and w' = 0 | exact to roundoff |
| Crude ring-of-force statics for the clamp moment | 2.2 percent |
| Linear scaling of stress and deflection with load | 1e-9 relative |
| Cubic scaling of deflection with thickness | 1e-6 relative |
| Region splitting does not change the answer | 1e-9 relative |
| Zero load gives identically zero result | exact |

The ring-of-force estimate treats the thrust as a line load at the band centroid,
radius 42.5 mm, giving M = P (r_c - b) / (2 pi b) = 43.01 N.m/m against the
solver value 42.08 N.m/m.

## Known limitations

- The 8 mm cable hole is not modelled. It is a local feature; a stress
  concentration factor should be applied separately rather than resolved here.
- Transverse shear deformation is neglected. With overhang/thickness = 3.3 this
  is a second-order effect at these stress levels but should be confirmed.
- The load is assumed uniform over the epoxy band. The real epoxy bond stress
  distribution is not known and is likely non-uniform.
- The bracket is treated as rigid. Real bracket compliance would reduce the
  clamp moment.
- Body forces and the plate self-weight are neglected, as agreed with the
  project supervisor.

## Files

| File | Purpose |
| :---- | :---- |
| `plate.py` | Main solver |
| `geometry.py` | Radii, areas and section properties |
| `material.py` | Al 6082-T6 properties |
| `loads.py` | Thrust values and the epoxy band load case |
| `beam.py` | Clamped-strip comparison baseline |
| `run.py` | Prints the summary and writes the stress plot |
| `test_plate.py` | Test suite |

## Running

    python3 run.py
    .venv/bin/python -m pytest test_plate.py -q

The virtual environment is created with `--system-site-packages` so that numpy
and sympy remain available:

    python3 -m venv --system-site-packages .venv
    .venv/bin/pip install pytest
