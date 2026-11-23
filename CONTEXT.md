# Context

## Objective

Analytical structural analysis of the flat annular thrust plate of a high-power rocket,
accurate to within a few percent of ANSYS, as the basis for a later geometry optimisation.
Performed as an independent analysis with a university FEA professor, to be included in
the CV.

## CV deliverables

"""
title: "Rocket Thrust-plate Analytical Optimisation",
  [Built a Python solver for stress and deflection of part from analytical elasticity theory, validated within [X]% of ANSYS],
  [Optimised the ring geometry to cut mass by [X]% ([Y] g to [Z] g) at a safety factor of [N] under a [thrust, N] load],
"""

## Meta project organization

- Perform a Git commit on every semi-important change.
- If you have any questions, or anything is unclear: ask the user before assuming
- Validated results with ANSYS will be provided later on, and you will see a
  ANSYS_RESULTS.txt file in the root directory with data to reference, until then, use
  whatever tools you have at your disposal to validate the code, but you won't be able to
  validate the physics until I provide that file.
- You should code the solver in Python unless you have a very good reason to use another
  language
- The code should be very simple, structured, very easily readable, and without any python
  magic or clever one-liners. You should not write any comments.
- The main solver should fit in one file, you may use other files for tests or other stuff.

## Status

The solver, methodology and results are documented in `README.md`. The technical model
description previously in this file has been removed to avoid duplicating the README,
which is now the single source of truth for the method and the numbers.

## Book references

The book, Martin H. S. *Elasticity*, is in the root directory. It is very long, so only the
appropriate sections should be read.

Sadd contains no plate bending theory. Its two-dimensional solutions are plane stress and
plane strain, which are in-plane problems, whereas this part is loaded transversely. The
plate formulation used in the solver is therefore not sourced from that text, and the
relevant overlap is limited to the following.

| Section | Use |
| :---- | :---- |
| 3.8 | Stress relations in cylindrical coordinates |
| 7.6 | Polar coordinate formulation, strain displacement and equilibrium |
| 8.4.1, Example 8.6 | Lame-type two-constant boundary condition structure |

## Questions for Mertol

- Ask him to explain how Gauss theorem is used, and why the volume integrals of energy and
  work are easier to compute than the stresses directly.
- Confirm that Sadd has no transverse plate bending content, and agree which reference the
  plate formulation should be attributed to.

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

- Awaiting ANSYS_RESULTS.txt to validate the physics before any result is accepted. The
  optimised thickness leaves a safety factor of 2.5, so there is much less margin for
  modelling error than the original 6 mm plate had.
- Decide whether to apply a stress concentration factor for the 8 mm cable hole.
- Model the epoxy clamp over 43 - 55 mm as a distributed elastic restraint rather than a
  rigid edge. The two rigid limits currently reported bracket the answer.
- Check transverse shear, since the free span over thickness ratio is low.
- Confirm whether the optimised thickness is manufacturable, and whether the bracket and
  epoxy interfaces assume a minimum bearing thickness.
- Buckling of the thin section has not been checked.
- The 441.3 g weight originally tabulated for the part is inconsistent with a flat 26 mm
  wide, 6 mm thick aluminium annulus, which weighs 111.2 g. It was removed from the tables
  but the discrepancy has never been explained.
