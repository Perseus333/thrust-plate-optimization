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
- Model the epoxy clamp over 43 - 55 mm as a distributed elastic restraint rather than a
  rigid edge. The two rigid limits currently reported bracket the answer.
- Check transverse shear, since the free span over thickness ratio is low.
- Buckling of the thin section has not been checked.
