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

# Open items

- Awaiting ANSYS_RESULTS.txt to validate the physics before any result is accepted. The
  optimised thickness leaves a safety factor of 2.5, so there is much less margin for
  modelling error than the original 6 mm plate had.
- The clamp is modelled as a single rigid edge at the middle radius of the epoxy band,
  r = 49 mm, which removed the ambiguity this item used to carry. If the epoxy band later
  turns out to be soft, the richer alternative is to clamp the plate at r = 43 mm and add
  a distributed elastic restraint over 43 - 55 mm; that needs an assumed epoxy stiffness.
  The two-clamp bracket showed the modelling choice moves the answer by tens of percent,
  so this stays the largest single uncertainty in the model.

# Verification performed

- The solution basis, the moment, shear and stress formulas, the continuity conditions at
  r = 36 mm and the boundary conditions at r = 29 mm and r = 49 mm were checked
  independently: residuals of the linear system are at machine precision, the clamp line
  reaction reproduces the applied load to twelve figures, and the clamped-edge hoop
  moment is exactly the Poisson ratio times the radial moment.
- The stress and deflection were reproduced by two independent routes, a from-scratch
  assembly of the same analytical solution and a Rayleigh-Ritz energy minimisation with an
  unrelated polynomial basis. They agree with the solver to 1e-13 and 0.12 percent
  respectively.
- Still outstanding: comparison against ANSYS once ANSYS_RESULTS.txt arrives. This has not
  been done and no result should be treated as validated until it is.

# Closing in this pass

- Clamp ambiguity over 43 - 55 mm: replaced by the single middle-radius clamp above.
- Transverse shear: checked at the ultimate load. The plate is moderately thick: the span
  from the load edge to the clamp is 13 mm, 2.17 times the 6 mm thickness and 4.38 times
  the optimised 2.97 mm, though the plain free band of 36 - 43 mm is only 7 mm wide. The
  shear share of the deflection falls as the plate thins and the peak shear stress is
  7.05 MPa against 100 MPa of bending stress. Transverse shear is not governing.
- Buckling: checked at the ultimate load. The critical pressure of the annulus is 29.6 MPa
  against an applied 2.21 MPa, a factor of 13.4, and the thin section is far too short to
  act as a column. The narrowest margin is the hoop compression under the load band,
  whose critical force 37.0 kN exceeds the 15.4 kN carried there by a factor of 2.4, the
  same order as the bending safety factor. Elastic instability is not governing either,
  but it is no longer unchecked.

