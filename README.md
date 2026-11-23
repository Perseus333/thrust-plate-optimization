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

## Optimisation

The plate must reach a safety factor of 2.5 against yield at the 2.5x ultimate
thrust of 3152 N, so the allowable peak stress is 250 / 2.5 = 100 MPa.

Because the peak bending stress of a plate scales as 1 / t^2, the required
thickness follows in closed form, without any search:

    t_required = t_reference * sqrt(sigma_reference / (yield / N))

Using the 6 mm plate, whose peak stress at 3152 N is 15.474 MPa:

    t_required = 6 * sqrt(15.474 / 100) = 2.3602 mm

| Quantity | 6 mm as built | 2.360 mm optimised |
| :---- | :---- | :---- |
| Peak von Mises at 3152 N | 15.47 MPa | 100.00 MPa |
| Safety factor at 3152 N | 16.16 | 2.50 |
| Peak von Mises at 1261 N | 6.19 MPa | 40.00 MPa |
| Safety factor at 1261 N | 40.39 | 6.25 |
| Maximum deflection at 3152 N | 5.25 um | 86.30 um |
| Annulus mass | 111.152 g | 43.724 g |

Mass saved is 67.428 g, or 60.66 percent. This is the mass of the flat annulus
alone, from the 110 mm outer diameter, 58 mm inner diameter and the stated
thickness. The 441.3 g weight originally tabulated for the part is not consistent
with such an annulus, which weighs 111.2 g at 6 mm, so it is not used as the
baseline.

The thinner plate is still a valid thin plate. Deflection falls from 0.09 percent
to 3.7 percent of the thickness, which remains far inside the small deflection
regime, and both overhang over thickness and radius over thickness improve. The
matrix conditioning also improves, from 2.4e6 to 1.4e5.

Stress at the optimised thickness is concentrated near the clamp and decays to
7.4 percent of peak at the free edge, so a tapered plate would save more. A
thickness profile scaled as 1 / sqrt(sigma) would take the plate to about 23.5 g,
a further 46 percent, but this needs the solver extended to variable thickness and
a large deflection analysis at the thin outer edge before it can be trusted.


| File | Purpose |
| :---- | :---- |
| `plate.py` | Main solver |
| `geometry.py` | Radii, areas and section properties |
| `material.py` | Al 6082-T6 properties |
| `loads.py` | Thrust values and the epoxy band load case |
| `optimize.py` | Thickness optimiser |
| `beam.py` | Clamped-strip comparison baseline |
| `run.py` | Prints the summary and writes the stress and mass plots |
| `test_plate.py` | Test suite |

## Running

    python3 run.py
    .venv/bin/python -m pytest test_plate.py -q

The virtual environment is created with `--system-site-packages` so that numpy
and sympy remain available:

    python3 -m venv --system-site-packages .venv
    .venv/bin/pip install pytest
