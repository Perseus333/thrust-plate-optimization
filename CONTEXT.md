# Context

## Objective:
Trying to perform analytical structural analysis of flat annular plate with accuracy comparable to ANSYS simulations to perform topological optimization for the given load.

## Context:
to be used as the thrust plate of a high-power rocket. The load is applied in the normal direction of the flat plate. To be simplified by simplifying through radial symmetry, converting it into a 1D beam. Clamped at one end through a bracket, and force applied at the other end through an expoxy region carrying the thrust from the rocket. 

Performed as independent analysis with FEA professor from university, to be included in the CV.

Deliverables for the CV entry:
"""
title: "Rocket Thrust-plate Analytical Optimisation",
  [Built a Python solver for stress and deflection of part from analytical elasticity theory, validated within [X]% of ANSYS],
  [Optimised the ring geometry to cut mass by [X]% ([Y] g to [Z] g) at a safety factor of [N] under a [thrust, N] load],
"""


## Tabulated data:

Thrust plate characteristics:

| Property​ | Value​ |
| :---- | :---- |
| Outer Diameter​ | 110 mm​ |
| Inner Diameter​ | 58 mm​ |
| Cable Hole Diameter​ | 8 mm​ |
| Bracket Outer Diameter​ | 70 mm​ |
| Thickness​ | 6 mm​ |
| Material​ | Al 6082 T6​ |
| Total Weight of the rocket​ | \~ 11 000 g |


Rocket motor characteristics:

| Property​ | Value​ |
| :---- | :---- |
| Peak Thrust​ | 1261 N​ |
| Peak Thrust, 2.5 safety​ | **3152 N**​ |
| Burn Time​ | 2.7 s​ |
| Total Impulse​ | 2714 Ns​ |


## Meta project organization

- Perform a Git commit on every semi-important change.
- If you have any questions, or anything is unclear: ask the user before assuming
- Validated results with ANSYS will be provided later on, and you will see a ANSYS_RESULTS.txt file in the root directory with data to reference, until then, use whatever tools you have at your disposal to validate the code, but you won't be able to validate the physics until I provide that file.
- You should code the solver in Python unless you have a very good reason to use another language
- The code should be very simple, structured, very easily readable, and without any python magic or clever one-liners. You should not write any comments.
- The main solver should fit in one file, you may use other files for tests or other stuff.

# Reading & Formulas to be used where necessary

The book can be found in the root directory of the workspace. It is very long, you should only read the appropriate sections if necessary.

## Elasticity by Martin H.S.

### Questions for Mertol
- Ask him to explain how Gauss theorem is used, and why the volume (energies/work) are easier/given
### Chapters to read
> [!expand]-
> (From [Gemini](https://gemini.google.com/app/d7d08db482ef9989))
> **Chapter 1: Mathematical Preliminaries**
> * **Must Read:**
> * **1.1, 1.2, 1.3:** Index notation and Kronecker delta ($\delta_{ij}$). Indicial notation is used throughout FEA literature and differential equations.
> * **1.8.1:** Divergence / Gauss Theorem. Crucial for converting strong-form differential equations into weak-form integral equations (the core of FEA stiffness matrices).
> 
> 
> * **Skim:**
> * **1.6:** Principal values (eigenvalues). Relevant later for post-processing, not initial FEA matrix setup.
> * **1.7:** Tensor algebra basics.
> 
> 
> * **Skip for now:** 1.4, 1.5, 1.8.2–1.8.4, 1.9.
> 
> **Chapter 2: Deformation: Displacements and Strains**
> 
> * **Must Read:**
> * **2.1, 2.2:** Small deformation theory and strain-displacement equations ($\epsilon = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T)$). Directly forms the basis for the FEA Strain-Displacement Matrix $[B]$.
> * **2.7:** Cylindrical coordinates ($r, \theta, z$). Required for ring geometry formulation.
> 
> 
> * **Skim:** 2.3, 2.4, 2.5.
> * **Skip for now:** 2.6 (Strain compatibility is guaranteed automatically when using continuous FEA shape functions).
> 
> **Chapter 3: Stress and Equilibrium**
> 
> * **Must Read:**
> * **3.1, 3.2:** Traction vectors and stress tensor ($\sigma_{ij}$).
> * **3.7:** Equilibrium equations ($\sigma_{ij,j} + F_i = 0$). Represents the strong-form governing equations solved by FEA.
> * **3.8:** Stress relations in cylindrical coordinates ($r, \theta, z$).
> 
> 
> * **Skim:** 3.5 (Von Mises stress, useful later for calculating safety factors).
> * **Skip for now:** 3.3, 3.4, 3.6.
> 
> **Chapter 4: Material Behavior—Linear Elastic Solids**
> 
> * **Must Read:**
> * **4.1, 4.2:** Generalized Hooke’s Law ($\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\epsilon}$). Defines the Constitutive Matrix $[D]$ used in every FEA solver.
> * **4.3:** Physical meaning of Elastic Modulus ($E$), Poisson's Ratio ($\nu$), Shear Modulus ($G$).
> 
> 
> * **Skip for now:** 4.4 (Thermoelasticity, unless the ring experiences thermal loads).
> 
> **Chapter 5: Formulation and Solution Strategies**
> 
> * **Must Read:**
> * **5.1:** Summary of field equations (brings Ch. 2, 3, and 4 together).
> * **5.2:** Boundary conditions (essential displacement BCs vs. natural traction BCs). Setting these wrong causes singular matrix errors in code.
> * **5.4:** Displacement formulation. Standard FEA is displacement-based.
> * **5.7.5, 5.7.6:** Approximate and numerical solution procedures (bridges elasticity theory to Rayleigh-Ritz/Galerkin concepts).
> 
> 
> * **Skim:** 5.3, 5.5, 5.6.
> * **Skip for now:** 5.7.1–5.7.4 (Analytical, inverse, and semi-inverse manual solution methods).

## Meeting 2026-09-09
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
