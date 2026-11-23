# Context

## Objective

Analytical structural analysis of the flat annular thrust plate of a high-power
rocket, accurate to within a few percent of ANSYS, as the basis for a later
geometry optimisation. Performed as an independent analysis with a university
FEA professor, to be included in the CV.

## CV deliverables

"""
title: "Rocket Thrust-plate Analytical Optimisation",
  [Built a Python solver for stress and deflection of part from analytical elasticity theory, validated within [X]% of ANSYS],
  [Optimised the ring geometry to cut mass by [X]% ([Y] g to [Z] g) at a safety factor of [N] under a [thrust, N] load],
"""

## Meta project organization

- Perform a Git commit on every semi-important change.
- If you have any questions, or anything is unclear: ask the user before assuming
- Validated results with ANSYS will be provided later on, and you will see a ANSYS_RESULTS.txt file in the root directory with data to reference, until then, use whatever tools you have at your disposal to validate the code, but you won't be able to validate the physics until I provide that file.
- You should code the solver in Python unless you have a very good reason to use another language
- The code should be very simple, structured, very easily readable, and without any python magic or clever one-liners. You should not write any comments.
- The main solver should fit in one file, you may use other files for tests or other stuff.

# Technical definition

## The part

A flat annular plate, clamped around its inner region by a bracket and loaded in
the normal direction by thrust carried through an epoxy band.

| Property | Value |
| :---- | :---- |
| Outer Diameter | 110 mm |
| Inner Diameter | 58 mm |
| Cable Hole Diameter | 8 mm |
| Bracket Outer Diameter | 70 mm |
| Thickness | 6 mm |
| Material | Al 6082 T6 |
| Total Weight of the rocket | ~ 11 000 g |

Rocket motor characteristics:

| Property | Value |
| :---- | :---- |
| Peak Thrust | 1261 N |
| Peak Thrust, 2.5 safety | 3152 N |
| Burn Time | 2.7 s |
| Total Impulse | 2714 Ns |

## Derived geometry

| Quantity | Value |
| :---- | :---- |
| Ring radial width, (110 - 58) / 2 | 26 mm |
| Clamp radius, bracket OD / 2 | 35 mm |
| Free outer edge radius | 55 mm |
| Unsupported overhang | 20 mm |
| Load band (epoxy) | 35 to 50 mm |
| Bearing area | 4005.5 mm^2 |
| Equivalent pressure at 1261 N | 314.8 kPa |
| Ratio of outer radius to thickness | 9.17 |

## Model

The part is a closed annular plate in bending, not a cantilever beam. Radial
symmetry reduces the problem to one coordinate, r, but the load path runs around
the circumference in hoop bending, which a radial strip cannot represent.

Three points distinguish this from the earlier 1D beam simplification:

1. The clamp acts at r = 35 mm, the bracket edge, not at the inner bore r = 29 mm.
2. The 26 mm figure in the tables is the radial width of the ring, not a beam
   length. The 5 mm and 15 mm split in the original drawing measured from the
   bracket edge, so the unsupported overhang is 20 mm.
3. The peak bending moment occurs at the clamped edge, which is exactly where a
   simple beam analogy is least reliable.

The clamped-strip baseline is retained in `beam.py` for comparison only. It
over-predicts the peak stress by a factor of about 26, and reports a safety factor
of 1.55 where the plate gives 40.4.

## Governing equations

Axisymmetric Kirchhoff thin plate bending. The deflection satisfies

    D * laplacian(laplacian(w)) = q(r),      D = E t^3 / (12 (1 - nu^2))

The homogeneous solution is

    w(r) = C1 + C2 r^2 + C3 log(r) + C4 r^2 log(r)

and the particular integral for uniform band pressure is

    w_p(r) = q r^4 / (64 D)

Moments and shear follow from the curvature relations

    M_r = -D ( w'' + nu w' / r )
    M_t = -D ( w' / r + nu w'' )
    Q   = -D d/dr ( w'' + w' / r )

Bending stresses are recovered from

    sigma_r = 6 M_r / t^2
    sigma_t = 6 M_t / t^2

The domain is split into regions at every load band discontinuity. Each region
carries the homogeneous form plus its own particular integral, and regions are
joined by continuity of w, w', M_r and Q. The four remaining conditions are

- clamped edge, r = 35 mm:  w = 0 and w' = 0
- free edge, r = 55 mm:     M_r = 0 and Q = 0

This yields a small dense linear system which is solved directly. The derivative
needed for the shear condition is available in closed form, since for the
homogeneous basis

    d/dr ( w'' + w' / r ) = 4 C4 / r

so no numerical differentiation is used anywhere in the solver.

## Results at 1261 N

| Quantity | Value |
| :---- | :---- |
| Peak von Mises stress | 6.19 MPa at r = 35 mm |
| Peak radial stress | -7.01 MPa |
| Peak hoop stress | -2.31 MPa |
| Maximum deflection | 2.10 um at r = 55 mm |
| Safety factor against yield, 250 MPa | 40.4 |
| Safety factor against ultimate, 300 MPa | 48.5 |

At the 2.5x load of 3152 N the peak stress is 15.47 MPa, a safety factor of 16.2
against yield. Peak deflection is 0.035 percent of the plate thickness, so the
small deflection and linearity assumptions hold with a very large margin.

## Validation status

ANSYS results are not yet available, so the solver is validated against
independent analytical checks. All 27 tests in `test_plate.py` pass.

| Check | Agreement |
| :---- | :---- |
| Exact symbolic annulus solution, uniform load | 1e-9 relative |
| Global equilibrium, clamp reaction equals applied load | 1e-5 relative |
| Free edge conditions M_r = 0 and Q = 0 | exact to roundoff |
| Clamped edge conditions w = 0 and w' = 0 | exact to roundoff |
| Crude ring-of-force statics for the clamp moment | 2.2 percent |
| Linear scaling of stress and deflection with load | 1e-9 relative |
| Cubic scaling of deflection with thickness | 1e-6 relative |
| Region splitting does not change the answer | 1e-9 relative |

## Known limitations

- The 8 mm cable hole is not modelled. It is a local feature, so a stress
  concentration factor should be applied separately.
- Transverse shear deformation is neglected. With overhang over thickness of 3.3
  this is a second-order effect at these stress levels, but should be confirmed.
- The load is assumed uniform over the epoxy band; the real bond stress
  distribution is likely non-uniform.
- The bracket is treated as rigid. Real compliance would reduce the clamp moment.
- Body forces and plate self-weight are neglected, as agreed with the supervisor.

# Optimisation findings

## Target

The 2.5 safety factor is applied to the peak thrust, giving an ultimate load of
3152 N, and the plate must reach a safety factor of 2.5 against yield at that
load. The allowable peak stress is therefore

    250 MPa / 2.5 = 100 MPa

## Scaling law

The peak bending stress of a plate is proportional to the bending moment divided
by the square of the thickness, so at fixed load

    sigma ~ 1 / t^2

and therefore, for a required safety factor N,

    t ~ sqrt(N)

The required thickness follows directly without any search:

    t_required = t_reference * sqrt(sigma_reference / (yield / N))

Using the 6 mm plate, whose peak stress at 3152 N is 15.474 MPa, gives

    t_required = 6 * sqrt(15.474 / 100) = 2.3602 mm

This closed form is verified to a safety factor of 2.5000000 against the solver.

## Result

| Quantity | 6 mm as built | 2.36 mm optimised |
| :---- | :---- | :---- |
| Thickness | 6.000 mm | 2.360 mm |
| Peak von Mises at 3152 N | 15.47 MPa | 100.00 MPa |
| Safety factor at 3152 N | 16.16 | 2.50 |
| Peak von Mises at 1261 N | 6.19 MPa | 40.00 MPa |
| Safety factor at 1261 N | 40.39 | 6.25 |
| Maximum deflection at 3152 N | 5.25 um | 86.30 um |
| Annulus mass | 111.152 g | 43.724 g |

    mass saved = 67.428 g
    mass saving = 60.66 percent

The mass figures are the mass of the annular plate alone, from the outer diameter
of 110 mm, the inner diameter of 58 mm and the quoted thickness. The 441.3 g
weight originally listed for the part is not consistent with a flat 26 mm wide
6 mm thick aluminium annulus, which weighs 111.2 g, and was removed from the
tables by the user. It is therefore not used as the optimisation baseline, and
the saving quoted above is relative to the 6 mm flat plate.

## Validity of the optimised design

Making the plate thinner pushes the thin plate assumptions harder, so the
dimensionless ratios are rechecked at 2.360 mm.

| Ratio | 6 mm | 2.360 mm | Requirement |
| :---- | :---- | :---- | :---- |
| Deflection over thickness | 0.00088 | 0.0366 | much less than 1 |
| Overhang over thickness | 3.33 | 8.47 | greater than about 5 |
| Outer radius over thickness | 9.17 | 23.30 | greater than about 5 |

Both directional ratios improve, because a thinner plate is relatively wider.
The deflection ratio rises from 0.09 percent to 3.7 percent of the thickness but
remains far inside the small deflection regime, so the linear Kirchhoff model
remains valid. The matrix conditioning also improves, falling from 2.4e6 at 6 mm
to 1.4e5, so there is no numerical penalty.

The design is therefore valid under the current model, though it now has much
less margin against any effect the model omits. Transverse shear deformation,
a non-uniform epoxy bond stress, and bracket compliance would each reduce the
available margin, and at a safety factor of 2.5 there is far less room for them
than at 16. This is the main risk in accepting the 2.36 mm result before ANSYS
validation.

## Stress distribution and the case for a taper

At 2.360 mm and 3152 N the stress is strongly non-uniform across the width.

| Radius, mm | von Mises, MPa | Fraction of peak |
| :---- | :---- | :---- |
| 35.0 | 100.00 | 100.0 percent |
| 38.0 | 59.06 | 59.1 percent |
| 40.0 | 39.04 | 39.0 percent |
| 42.0 | 24.16 | 24.2 percent |
| 45.0 | 11.34 | 11.3 percent |
| 48.0 | 8.85 | 8.8 percent |
| 50.0 | 8.60 | 8.6 percent |
| 55.0 | 7.41 | 7.4 percent |

Only the material near the clamp is working hard. A uniform plate is therefore
wasteful, and a plate tapered in thickness would save considerably more mass.
Scaling the local thickness as 1 / sqrt(sigma) gives a required profile falling
from 2.360 mm at the clamp to 0.642 mm at the free edge, which would reduce the
mass of the 35 to 55 mm region from 36.0 g to about 15.8 g, a further 56 percent
of that region.

Over the whole plate the saving is smaller, because the region from 29 to 35 mm
lies under the bracket and must stay at the full 2.360 mm to carry the clamp
bearing. Including that 7.7 g, the tapered plate would weigh about 23.5 g against
43.7 g for the uniform plate, a further 46 percent.

That estimate is indicative only, because the solver currently assumes a constant
flexural rigidity, so it cannot yet model a taper. It is also optimistic: at
0.64 mm the local deflection reaches about 6.8 times the local thickness, which
breaks the small deflection assumption and would require a large deflection
analysis with membrane stiffening. A real tapered design would therefore land
somewhere between 23.5 g and 43.7 g, and needs the solver extended to variable
thickness before the figure can be trusted.

## Limitations of this optimisation

- Only the uniform thickness is optimised. The radii, the load band and the
  bracket diameter are unchanged.
- Only through thickness stress is considered. Deflection is not constrained, and
  at 86 um it is large in relative terms even though it is structurally harmless.
- Buckling of the thin section under thrust is not checked. With a 20 mm span and
  a 2.36 mm thickness this is not expected to govern, but it has not been verified.
- The cable hole stress concentration is not included.

# Book references

The book, Martin H. S. *Elasticity*, is in the root directory. It is very long,
so only the appropriate sections should be read.

Sadd contains no plate bending theory. Its two-dimensional solutions are plane
stress and plane strain, which are in-plane problems, whereas this part is loaded
transversely. The plate formulation above is therefore not sourced from that
text, and the relevant overlap is limited to the following.

| Section | Use |
| :---- | :---- |
| 3.8 | Stress relations in cylindrical coordinates |
| 7.6 | Polar coordinate formulation, strain displacement and equilibrium |
| 8.4.1, Example 8.6 | Lame-type two-constant boundary condition structure |

## Chapters of interest

> [!expand]-
> **Chapter 1: Mathematical Preliminaries**
> * **1.1, 1.2, 1.3:** Index notation and Kronecker delta.
> * **1.8.1:** Divergence / Gauss theorem, converting strong form to weak form.
> * **Skim:** 1.6, 1.7. **Skip:** 1.4, 1.5, 1.8.2 to 1.8.4, 1.9.
>
> **Chapter 2: Deformation: Displacements and Strains**
> * **2.1, 2.2:** Small deformation theory and strain displacement relations.
> * **2.7:** Cylindrical coordinates.
> * **Skim:** 2.3, 2.4, 2.5. **Skip:** 2.6.
>
> **Chapter 3: Stress and Equilibrium**
> * **3.1, 3.2:** Traction vectors and the stress tensor.
> * **3.7:** Equilibrium equations.
> * **3.8:** Stress relations in cylindrical coordinates.
> * **Skim:** 3.5. **Skip:** 3.3, 3.4, 3.6.
>
> **Chapter 4: Material Behavior, Linear Elastic Solids**
> * **4.1, 4.2:** Generalised Hooke's law.
> * **4.3:** Elastic modulus, Poisson's ratio, shear modulus.
> * **Skip:** 4.4.
>
> **Chapter 5: Formulation and Solution Strategies**
> * **5.1:** Summary of field equations.
> * **5.2:** Essential displacement and natural traction boundary conditions.
> * **5.4:** Displacement formulation.
> * **5.7.5, 5.7.6:** Approximate and numerical solution procedures.
> * **Skim:** 5.3, 5.5, 5.6. **Skip:** 5.7.1 to 5.7.4.

## Questions for Mertol

- Ask him to explain how Gauss theorem is used, and why the volume integrals of
  energy and work are easier to compute than the stresses directly.
- Confirm that Sadd has no transverse plate bending content, and agree which
  reference the plate formulation should be attributed to.

# Meeting 2026-09-09

> [!Notes meeting 2026-09-09]-
> questions for Mertol:
> - what is that makes Gauss theorem (aka the energy calculation) easy to compute?
> 	total poential energy and total kinetic energy - stress is the energy - first you calculate the energy and then you calculate the stress from that
> 	calculate derivative, find 0, find a local maximum/minimum - mechanical systems derivative is the energy - if slope is 0 then energy is minimized - most of the time it's the potential energy
> - do we usually assume small deformations?
> 	YES - if you can express the equilibrium on the deformed shape like the undeformed then it's small - no displacements
> 	if you assume deformedne shape and then express equilibrium that's very hard to solve
> - and when solving the equations why can't we use the engineering strain?
> 	handles mathematical problem mathematically, engineering stress is a simplification, it's empirical
> - what is the form that is most commonly used? derivatives, index notation, tensor notation? which one will be most useful to us?
> 	tensor notation - pytorch
> - is it true that the thrust ring will be multiply connected?
> 	whether it transfer the moment - multiply connected
> - at what ratio of diameter / thickness can we assume the structure to be 2D? is it 1/10 like in the beam examples from class?
> 	yes
> - I assume we would use the cylindrical coordinates?
> 	useful for arches, half space (inifite half 3D space) - business/bussiness problem 
> - is it safe to ignore body forces given the small weight of the component compared to the thrust of the motor?
> 	probably yes - acceleration
> - we would have mixed boundaries beause we have traction and displacement, no?
> 	mixed, because it can move and can have applied forces
> 	- if so, what are the two orthogonal directions that we can choose to analyze? x, y?
> 
> - can we assume fixed boundaries? or if they deform they cannot be fixed?
> 	yeah, 
> - which of the analytical procedures described in 5.7.4 should we use?
> 	we don't know yet
> - plane stress is what we should use, right?
> 	exactly
> - do we need to use navier equations? or simple equilibrium equations are fine enough?
> 	heavily mathematical, equilibrium equations are more physical based - use equilibrium equations; sometimes navier is better, and sometimes they are they only 
> - do we need to use polar coordinates if we simplify it down to 1D beams?
> 	probably not
> 
> 
> summarize my problem first, try to match with examples in the book
> domain
> equilibrium equations
> formulate the equilibrium
> unit equilibrium equation
> finding the governing equation
> 
> find similar examples in the book - overlapping, not necessarily the same
> 	try to mimic the steps in the book
> 
> 
> Research in the future: computationally efficient  	 galerkin method

# Open items

- Awaiting ANSYS_RESULTS.txt to validate the physics before the 2.360 mm result is
  accepted. At a safety factor of 2.5 there is much less margin for any modelling
  error than the original 6 mm plate had.
- Decide whether to apply a stress concentration factor for the 8 mm cable hole.
- Uniform thickness optimisation is complete. Extending the solver to variable
  thickness would allow a taper, estimated to save a further 56 percent, but
  requires a large deflection analysis at the thin outer edge.
- Confirm whether the 2.36 mm plate is manufacturable, and whether the bracket
  and epoxy interfaces assume a minimum bearing thickness.
- Buckling of the thin section has not been checked.
