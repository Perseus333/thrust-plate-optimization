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

- Awaiting ANSYS_RESULTS.txt to validate the physics.
- Decide whether to apply a stress concentration factor for the 8 mm cable hole.
- Optimisation of the ring geometry is not started. At a safety factor of 40 the
  current 6 mm plate is substantially over-designed, so there is significant mass
  to remove.
