# Rocket Thrust Plate — Analytical Stress, Deflection and Mass Optimisation

**Independent findings report for human review**

| | |
| :---- | :---- |
| Part | Flat annular aluminium thrust plate, high-power rocket |
| Analysis type | Axisymmetric thin-plate bending, closed-form semi-analytical |
| Tooling | Python 3.14, NumPy, Matplotlib, pytest |
| Status | Solver complete and self-consistent; **physics not yet validated against ANSYS** |
| Report date | 2026-11-23 |

---

## 1. What this report is

This document explains, in plain language, what was analysed, how it was analysed, what
was found, and what is still uncertain. It is written so that a reader with general
engineering literacy — but no prior exposure to this project or to plate theory — can
follow the whole argument.

Every number quoted here was **re-derived from the code during the preparation of this
report**, not copied from notes. Where the project's own notes disagree with what the
code actually produces, the disagreement is flagged explicitly in
[Section 12](#12-discrepancies-and-issues-found-in-the-written-record).

The central caution, stated up front and repeated throughout:

> **The solver has never been checked against a real structural simulation.** It has been
> checked extensively against *itself* and against *hand calculations*, and it passes all
> of those checks. But those checks confirm that the code correctly solves the equations
> it was given. They do **not** confirm that those equations describe the real part. That
> question is still open, and it is the single most important limitation in this report.

---

## 2. Plain-language summary of the findings

The part is a flat ring of aluminium, 110 mm across, with a 58 mm hole in the middle, 6 mm
thick. A rocket motor pushes on it through a ring-shaped band of glue (epoxy), and a
bracket holds it near the middle. The questions were: how much stress does it see, how
much does it bend, and how much lighter could it be?

**Headline results:**

1. **The existing 6 mm plate is extremely lightly loaded.** At the motor's peak thrust of
   1261 N, the highest stress anywhere in it is about **6.19 MPa**. The aluminium it is
   made from (6082-T6) starts to yield permanently at roughly **250 MPa**. That is a
   safety factor of about **40**. The plate barely notices the load.

2. **It barely bends.** The maximum movement is **2.10 micrometres** — about one fiftieth
   of the thickness of a human hair, and roughly 0.035% of the plate thickness. It is,
   for all practical purposes, perfectly rigid.

3. **It is therefore much heavier than it needs to be.** To meet a target safety factor
   of **2.5** at 2.5× the peak thrust (3152 N), the plate only needs to be
   **2.36 mm thick** instead of 6 mm.

4. **That saves 60.7% of the mass** — from **111.15 g down to 43.72 g**. A saving of
   **67.43 g** on a rocket that weighs about 11 kg. This is a small absolute saving
   (about 0.6% of the rocket's mass) but a large *relative* one on this component.

5. **The stress is very unevenly distributed.** At the 2.36 mm thickness, the stress is
   100 MPa right at the clamped edge and falls to 7.4 MPa at the outer rim. The outer
   two-thirds of the plate is doing almost no structural work. This strongly suggests a
   *tapered* plate could save even more weight — but the current solver cannot model a
   taper, so any number for that is an estimate only.

6. **The 2.36 mm result carries real risk.** At a safety factor of 2.5, there is very
   little room left for any effect the model ignores. Several known omissions (described
   in [Section 10](#10-limitations)) would each *reduce* the true safety margin. Until
   ANSYS confirms the physics, **2.36 mm should be treated as a promising candidate, not
   an accepted design.**

A parallel finding worth noting: a first attempt that treated the part as a simple
cantilever beam over-predicted the stress by a factor of about **26**. That is discussed
in [Section 6](#6-progression-of-methods), because understanding *why* it was wrong is
central to trusting the current model.

---

## 3. The physical problem

### 3.1 The part

A flat annular plate — an "annulus", meaning a ring, i.e. a disc with a concentric hole in
it. Think of a washer or a donut, but flat and thin.

| Property | Value |
| :---- | :---- |
| Outer diameter | 110 mm |
| Inner diameter (central bore) | 58 mm |
| Cable hole diameter | 8 mm |
| Bracket outer diameter | 70 mm |
| Thickness | 6 mm |
| Material | Aluminium 6082-T6 |
| Total rocket mass | ~11 000 g |

The **cable hole** is a small 8 mm hole somewhere in the ring, presumably for routing a
wire. It is not modelled — see [Section 10](#10-limitations).

### 3.2 How the load arrives and how it is held

Two things act on the plate:

- **The motor pushes.** The rocket motor's thrust is transferred into the plate through a
  ring-shaped band of **epoxy** (structural adhesive). That band spans from radius 35 mm
  to radius 50 mm. Because the band is an annulus of known area, the total force can be
  converted into an equivalent **pressure** — force per unit area — spread over the band.
- **The bracket holds.** A bracket clamps the plate around its inner region. The bracket's
  outer diameter is 70 mm, so it grips out to a radius of **35 mm**.

### 3.3 The single most important geometric distinction

This caused the greatest confusion in the early work, so it is worth stating very
carefully.

| Quantity | Value | What it actually is |
| :---- | :---- | :---- |
| Radial width of the ring | 26 mm | (110 − 58) / 2 — the *material width*, from bore to rim |
| Clamp radius | 35 mm | Where the bracket's grip ends |
| Free outer edge radius | 55 mm | The rim of the plate |
| **Unsupported overhang** | **20 mm** | 55 − 35 — the part that is genuinely free to bend |
| Load band (epoxy) | 35 to 50 mm | Where the motor's force enters |

The trap: **26 mm looks like a beam length, but it is not.** It is the width of material
from the central bore to the outer rim. Because the bracket grips all the way out to
radius 35 mm, the bore-to-bracket region (radii 29–35 mm) is clamped and cannot move. The
portion that can actually flex is only the **20 mm** from radius 35 mm to radius 55 mm.

An early version of the project used 26 mm as the bending span. It should have been 20 mm,
measured from the bracket edge. Separately, the original drawing's "5 mm and 15 mm" split
was also measured from the bracket edge, which is what establishes the 20 mm overhang.

Additionally — and this is the conceptual crux of the whole project — **the plate is not a
beam at all.** It is a closed ring. When it bends, the load path runs *around* the
circumference as well as radially. A narrow radial strip cut out of the ring cannot
represent that circumferential (hoop) action. This is exactly why the beam model failed so
badly.

### 3.4 The loads

| Property | Value |
| :---- | :---- |
| Peak thrust | 1261 N |
| Peak thrust × 2.5 safety | 3152 N |
| Burn time | 2.7 s |
| Total impulse | 2714 N·s |
| Bearing area (epoxy band) | 4005.53 mm² |
| Equivalent pressure at 1261 N | 314.81 kPa |

Derived quantities confirmed by direct computation:

| Quantity | Value |
| :---- | :---- |
| Annulus area (the ring material itself) | 6861.24 mm² |
| Mass of 6 mm annulus | 111.1521 g |
| Mass of 2.3602 mm annulus | 43.7235 g |
| Ring width / overhang | 26.00 mm / 20.00 mm |

The distinction between **bearing area** (4005.53 mm², the epoxy band) and **annulus area**
(6861.24 mm², the aluminium ring) matters and is a common source of confusion. They are
different regions and must not be interchanged.

### 3.5 Two load cases, and why both matter

- **1261 N** — the real peak thrust. This tells you how the part behaves in service.
- **3152 N** — 1261 N multiplied by 2.5. This is not a load that will ever occur in
  flight. It is a **design load**: the requirement is that the plate survives 2.5× the
  worst expected load without yielding. Sizing against 3152 N builds in the safety margin.

A subtlety that is easy to get wrong: the "2.5 safety factor" is applied **twice over**, in
two different roles. It is used once to *create* the 3152 N design load, and then the
design target is to achieve a safety factor of 2.5 *at that already-inflated load*. In
practical terms the requirement reduces to: the peak stress at 3152 N must not exceed
250 / 2.5 = **100 MPa**. That is the number the optimisation drives towards.

---

## 4. Jargon and notation explained

This section defines every technical term used elsewhere in the report. Nothing here
requires prior knowledge.

### 4.1 Structural mechanics terms

**Stress (σ, sigma).** Internal force per unit area inside a material, measured in pascals
(Pa). 1 MPa = 10⁶ Pa. A stress of 250 MPa means each square millimetre of material carries
250 newtons of internal load.

**Bending stress.** When a flat plate is pushed sideways (perpendicular to its face), it
curves. The outer surface on the convex side stretches, the other side compresses. This
stretching/compressing is bending stress. It is largest at the surfaces and zero at the
mid-plane.

**Deflection (w).** How far a point on the plate moves out of its original plane. Measured
here in micrometres (µm); 1 µm = 0.001 mm.

**Yield stress.** The stress above which a metal stops springing back and starts deforming
permanently. For 6082-T6 aluminium this is taken as **250 MPa**.

**Ultimate stress.** The stress at which the material actually breaks. Taken as **300 MPa**
here. Less relevant, because a part that has yielded is already a failure for a rocket
component.

**Safety factor.** Yield stress divided by the actual peak stress. A safety factor of 40
means the part could carry 40 times the current load before yielding. Higher is safer;
higher is also heavier, because extra material is usually what buys the margin.

**Poisson's ratio (ν, nu).** When you stretch a material in one direction it usually
contracts slightly in the perpendicular direction. Poisson's ratio is the ratio of that
sideways contraction to the applied stretch. For aluminium, ν ≈ **0.33**. It matters in
plate bending because it couples the radial and circumferential bending.

**Young's modulus (E).** A measure of stiffness — how much stress is needed to produce a
given stretch. For aluminium, E = **70 GPa** (70 × 10⁹ Pa). A higher E means a stiffer
material. It governs deflection strongly but affects *stress* barely at all in this
problem, because the load is fixed and the geometry is fixed.

**Density.** Mass per unit volume. Aluminium: **2700 kg/m³**.

### 4.2 Plate theory terms

**Annulus / annular plate.** A flat ring: a disc with a concentric hole. "Annular" simply
means ring-shaped.

**Axisymmetric.** Symmetric about an axis. Here: the part, the load and the supports are
all circularly symmetric, so nothing varies with angle around the circle. This is a huge
simplification — it reduces a two-dimensional problem to a **one-dimensional** one, where
everything depends only on radius *r*. This is why the solver is fast and exact rather than
approximate.

**Kirchhoff thin-plate theory.** The standard classical theory for bending of thin plates.
Its central assumption is that a line drawn perpendicular to the plate's mid-surface stays
perpendicular to it after bending, and does not change length. This is the plate equivalent
of the "plane sections stay plane" assumption in beam theory. It ignores deformation caused
by *shear* (see below). It is accurate when the plate is thin compared to its other
dimensions.

**Flexural rigidity (D).** The plate's resistance to bending, defined as

    D = E t³ / (12 (1 − ν²))

where t is thickness. The t³ is critical: doubling the thickness makes the plate **eight
times** stiffer. This single fact drives the entire optimisation.

**Bending moment (M_r, M_t).** The internal turning effect per unit length of a cross
section. Subscript *r* means **radial** — acting on a cut that runs around the circle.
Subscript *t* means **tangential** or **hoop** — acting on a cut that runs radially outward.

**Hoop stress / hoop bending.** Stress acting in the circumferential direction — i.e.
"around the ring", like the bands of a barrel. This is the ingredient a beam model cannot
capture, and its neglect is precisely why the beam model was wrong.

**Radial stress.** Stress acting along the radius — i.e. outward from the centre.

**Transverse shear (Q).** The force trying to slide one cross section past the next, like
scissors. In thin plates it is small, but it is what carries the load from where it is
applied to where it is reacted.

**von Mises stress.** A single combined number summarising a complex multi-directional
stress state, used to compare against the single yield stress from a tensile test. For a
plate with only radial and hoop stress, it is

    σ_vm = sqrt(σ_r² − σ_r σ_t + σ_t²)

This is the quantity compared against the 250 MPa yield stress to get a safety factor. It
is *not* simply the larger of the two stresses.

**Clamped edge.** A support that prevents both movement and rotation — like a beam concreted
into a wall. Written as w = 0 and dw/dr = 0 (zero deflection, zero slope).

**Free edge.** An edge with nothing attached to it, so no moment and no shear can be
transferred across it — written as M_r = 0 and Q = 0.

**Boundary conditions.** The constraints and applied forces at the edges of the region being
solved. Getting these right is most of the modelling work; getting them wrong is the most
common source of silently wrong answers.

**Stress concentration.** A local spike in stress caused by a geometric feature that
disrupts the smooth flow of internal force — a hole, notch, or sharp corner. The 8 mm cable
hole would cause one. Crucially, such a spike is a *local* effect: it does not change the
overall load path, so it is normally handled by multiplying the smooth-theory stress by a
factor near the hole, rather than by modelling the hole inside the global solver.

**Large deflection / membrane stiffening.** The classical thin-plate theory assumes the
deflection is small enough that the plate's geometry is essentially unchanged. When
deflections become comparable to the thickness, the plate starts carrying load by
*stretching* in its own plane rather than by pure bending, which makes it effectively
stiffer and redistributes stress. That regime needs a different, nonlinear theory. This is
why "deflection over thickness" is tracked throughout.

### 4.3 Numerical and software terms

**Linear system / matrix solve.** The physics reduces to a set of simultaneous linear
equations, written as a matrix times a vector of unknowns equals a vector of knowns. A
standard linear-algebra solver finds the unknowns exactly (to numerical precision) in one
step — no iteration.

**Matrix condition number.** A measure of how numerically delicate a linear system is. A
small condition number is good. A very large one (above roughly 10¹² in double precision)
means the answer is being computed from nearly-cancelling large numbers and may lose
accuracy. It is a diagnostic for numerical reliability, not a physical quantity.

**Particular integral.** For a differential equation with a forcing term, a specific
solution that reproduces the forcing. The **homogeneous solution** is the general solution
of the equation with no forcing. Adding the two gives the general solution of the forced
problem. The code uses one homogeneous solution per region plus that region's particular
integral.

**Regression test.** An automated check that a specific known-good behaviour still holds.
The project has 36 of them. They are the safety net that allows the code to be modified
without silently breaking previously correct behaviour.

**Semi-analytical.** A method that is analytical (closed-form mathematics, exact
expressions) in form, but requires a small numerical solve to pin down constants. This
solver is semi-analytical: it uses exact closed-form basis functions and solves only a
small set of coefficients numerically. It is not a finite element mesh, and it does not
approximate the stress field with piecewise polynomials — so it suffers none of the mesh
refinement error that a finite element model would.

**ANSYS.** Commercial finite element software, widely regarded as the industry standard.
It will eventually be used to check the physics independently. It is a completely different
method (numerical, mesh-based, three-dimensional), so agreement between the two would be
strong evidence that both are right.

**Finite element analysis (FEA).** A numerical technique that cuts a complex shape into many
small simple pieces, solves the physics on each, and stitches the results together. Very
general and very powerful, but approximate — the answer depends on how fine the mesh is.
The approach here is the opposite: exploit the symmetry to get exact mathematics, and use
the finite element result later as an independent check.

---

## 5. Method used

### 5.1 The governing equation

For an axisymmetric thin plate under a distributed transverse pressure q(r), the deflection
w(r) satisfies the Kirchhoff plate equation:

    D · ∇⁴w = q(r)

where ∇⁴ is the biharmonic operator and D is the flexural rigidity. In polar coordinates
with no angular variation, this becomes a fourth-order ordinary differential equation in
the single variable r.

Because it is **fourth order**, the general solution contains **four** arbitrary constants.
Those four constants are fixed by four boundary conditions — two at the clamped edge, two
at the free edge. This is the mathematical reason the model needs exactly the four
conditions listed below, no more and no fewer.

### 5.2 The homogeneous solution and the particular integral

With no load (q = 0), the equation is solved by the four-term basis:

    w(r) = C₁ + C₂ r² + C₃ ln(r) + C₄ r² ln(r)

Each term has a clear meaning: C₁ is a rigid-body translation; C₂ r² is a rigid-body-like
parabolic bowl; the two logarithmic terms are the genuinely "plate-like" bending shapes
concentrated near the hole. Note that ln(r) terms appear, which is exactly why the model
breaks down at r = 0 — but this part has a bore at 29 mm, so r is never zero and the
logarithms are perfectly well behaved.

For a **uniform** pressure q over the loaded band, a particular integral is

    w_p(r) = q r⁴ / (64 D)

### 5.3 Recovering moments and shear

Once w(r) is known, curvatures follow by differentiation, and the internal moments follow
from the curvature relations:

    M_r = −D ( w'' + ν w'/r )
    M_t = −D ( w'/r + ν w'' )
    Q   = −D d/dr ( w'' + w'/r )

and the surface bending stresses from

    σ_r = 6 M_r / t²
    σ_t = 6 M_t / t²

The 6/t² factor is the same section-modulus relation used in beam theory: it converts an
internal moment per unit width into a surface stress.

### 5.4 Handling the load discontinuity — the region method

The load is zero from radius 35 mm inward of the band, non-zero between 35 and 50 mm, and
zero again from 50 mm out to 55 mm. Each change is a **discontinuity** in the forcing, and a
discontinuous load produces a kink in the solution.

The solver handles this by splitting the domain at every load discontinuity into
**regions**, and giving each region its own four constants and its own particular integral.
Adjacent regions are then stitched together by requiring **continuity of w, w', M_r and Q**
across each shared boundary. This is the physically correct set of matching conditions: the
plate cannot tear apart (w continuous), cannot form a crease (w' continuous), and cannot
have a jump in internal moment or shear (M_r and Q continuous) at a location where nothing
is applied.

For the actual load case this produces **2 regions** (confirmed by direct inspection), and
hence an 8 × 8 linear system.

### 5.5 The four boundary conditions

- **Clamped edge, r = 35 mm:** w = 0 and w' = 0. The plate is held flat and cannot rotate.
  *This is the most important modelling decision in the project:* the clamp is placed at the
  bracket edge radius of 35 mm, **not** at the inner bore of 29 mm. Getting this wrong
  changes the answer substantially, because it changes how much of the plate is free to
  bend.
- **Free edge, r = 55 mm:** M_r = 0 and Q = 0. Nothing is attached to the rim, so no
  moment or shear can be transmitted through it.

### 5.6 Two details that add confidence

**No numerical differentiation anywhere.** The shear condition requires d/dr(w'' + w'/r).
For the homogeneous basis this derivative is available in exact closed form:

    d/dr ( w'' + w'/r ) = 4 C₄ / r

so the shear is computed exactly rather than by finite differences. Numerical
differentiation is a notorious source of error, and the solver avoids it entirely. (Note:
the *tests* use finite differences to check the solver independently — which is the right
way round. The solver is exact; the tests are deliberately independent.)

**Exact reduction to a known case.** When the load is applied uniformly over a single
region spanning the whole plate, the problem has a classical textbook closed-form solution.
The solver reproduces that solution to 1 part in 10⁹. This is a strong check that the
basis functions, the derivative relations and the assembly are all correct.

### 5.7 How the code is organised

The project follows a deliberately simple, readable structure with no comments and no
clever constructs, as required by the project brief.

| File | Purpose |
| :---- | :---- |
| `plate.py` | The main solver: basis functions, region assembly, solution, sampling, reporting |
| `geometry.py` | Radii, areas, section properties |
| `material.py` | Al 6082-T6 properties, flexural rigidity |
| `loads.py` | Thrust values, epoxy band, pressure conversion |
| `optimize.py` | Thickness optimiser |
| `beam.py` | The clamped-strip baseline, retained for comparison only |
| `run.py` | Prints the summary and writes both plots |
| `test_plate.py` | The 36-test suite |

Only `plate.py`, `optimize.py`, `geometry.py`, `material.py` and `loads.py` participate in
producing results. `beam.py` produces no result used in the final answer — it exists purely
as a documented comparison.

### 5.8 How the optimisation was performed

This is a point of elegance worth highlighting. No search, iteration or numerical
optimisation was needed.

For a plate in bending, the peak stress scales as the bending moment divided by thickness
squared. Since the moment itself is set by the fixed load and fixed geometry, at constant
load:

    σ ∝ 1 / t²

Therefore, for a required safety factor N against a fixed yield stress:

    t ∝ √N

and the required thickness follows **in closed form** from a single reference solve:

    t_required = t_reference · √( σ_reference / (yield / N) )

Substituting the 6 mm plate's peak stress at 3152 N (15.4741 MPa) and the 100 MPa
allowable:

    t_required = 6 · √(15.4741 / 100) = 2.3602 mm

This is not an approximation or a fit — it is a direct consequence of the governing
equation, and it is verified against the solver to give a safety factor of
**2.500000000000008**, i.e. 2.500 to far better than one part in 10¹⁴. Solving one load
case and scaling analytically is far more robust and far more transparent than running an
optimiser, and it makes the result easy to audit by hand.

---

## 6. Progression of methods

The final model was not the first attempt. The path taken is instructive, because the
failure of the early approach is what revealed the true nature of the problem.

### 6.1 Stage 1 — Beam idealisation (superseded)

**The reasoning.** The plate is circularly symmetric, so it was argued that a radial slice
could be treated as a one-dimensional cantilever beam: clamped at the inner end, loaded at
the outer end. This is a natural and common first simplification. The project notes record
a discussion with the supervisor confirming that at a diameter-to-thickness ratio of about
10, a two-dimensional idealisation was acceptable — and even that polar coordinates might
not be needed if the part were reduced to a beam.

**What it predicted** (at 1261 N, span 26 mm, width 26 mm):

| Quantity | Beam prediction |
| :---- | :---- |
| Root bending stress | 161.67 MPa |
| Tip deflection | 102.65 µm |
| Safety factor | 1.55 |

**Why it was wrong.** Three separate errors compounded:

1. **Wrong span.** The 26 mm ring width was used as the bending span. The real span is the
  20 mm unsupported overhang measured from the bracket edge.
2. **Wrong structural idealisation.** This is the fundamental error. A ring is not a strip.
   When an annulus bends, the material also bends *around* the circumference in hoop
   bending, and the ring derives a great deal of its stiffness from that. A radial strip
   cut from the ring has no hoop continuity and therefore no hoop stiffness — it is
   artificially floppy, so it over-predicts both stress and deflection.
3. **Wrong location for the critical stress.** The beam model puts the peak moment at the
   root, which happens to coincide in location with the plate's peak — but for entirely
   different reasons, and with a wildly different magnitude. The plate's peak moment arises
   from the restrained rotation at the clamp under axisymmetric bending, which a beam
   analogy cannot represent.

**How far off it was.** The beam model's 161.67 MPa against the plate model's 6.19 MPa is a
factor of **26.1**. The beam model reports a safety factor of 1.55 where the correct answer
is about 40. A design based on the beam model would have rejected a perfectly safe part,
or — worse, if used to size a plate — driven a needlessly heavy design.

This model is retained in `beam.py` purely as a documented cautionary baseline. It
contributes nothing to the final result.

### 6.2 Stage 2 — Axisymmetric plate formulation (current)

The model was rebuilt on the correct physical basis: axisymmetric Kirchhoff thin-plate
bending, solved semi-analytically by the region method described in
[Section 5](#5-method-used). The three corrections relative to the beam model were:

1. The clamp was moved to r = 35 mm, the bracket edge — not the inner bore at 29 mm.
2. The span was recognised as the 20 mm overhang, not the 26 mm ring width.
3. The circumferential (hoop) load path was included by solving the full axisymmetric plate
   equation rather than a radial strip.

This produced the final results: 6.19 MPa peak stress, 2.10 µm deflection, safety factor
40.4.

### 6.3 Stage 3 — Independent validation

Before trusting the new model, it was checked against a series of independent calculations.
These are catalogue in [Section 8](#8-validation-of-the-solver). Notably, the exploitation
of an exact symbolic annulus solution and the crude ring-of-force static estimate are
genuinely independent of the solver's own machinery.

### 6.4 Stage 4 — Optimisation

With a validated solver, the thickness was optimised in closed form (Section 5.8), giving
2.3602 mm and a 60.7% mass saving. The thin-plate validity ratios were then **re-checked at
the new thickness**, because making a plate thinner pushes the theory's assumptions harder
— a step that is easy to skip and important not to.

### 6.5 Stage 5 — Identification of the taper opportunity

Plotting the stress distribution at the optimised thickness revealed that stress decays
from 100 MPa at the clamp to 7.4 MPa at the rim. Most of the material is barely working.
This suggested a tapered plate, estimated (with strong caveats) at a further ~46% saving.
This stage is **analysis only** — the solver cannot yet model a variable thickness.

---

## 7. Findings: results at the peak thrust of 1261 N

| Quantity | Value |
| :---- | :---- |
| Peak von Mises stress | **6.1896 MPa** at r = 35.0 mm |
| Peak radial stress | **−7.0133 MPa** at r = 35.0 mm |
| Peak hoop stress | **−2.3144 MPa** at r = 35.0 mm |
| Maximum deflection | **2.1012 µm** at r = 55.0 mm |
| Safety factor against yield (250 MPa) | **40.39** |
| Safety factor against ultimate (300 MPa) | **48.47** |
| Free-edge hoop stress | −0.4586 MPa |
| Free-edge radial stress | 0.0000 MPa |

**Reading these numbers.**

- The peak stress occurs **exactly at the clamped edge, r = 35 mm**, and the radial and
  hoop peaks occur at the same place. This is the expected and physically sensible result
  for a clamped plate: the clamp restrains rotation, which generates the largest bending
  moment right at the restraint.
- The stresses are **negative**, meaning the loaded surface is in **compression**. For a
  plate pushed on one face while clamped at the inner edge, the surface being pushed is
  compressed near the fixed edge. The *sign* is a matter of sign convention and of which
  surface is examined; the **magnitude** is what matters for yield.
- The **radial stress dominates** (|7.01| vs |2.31| MPa). The hoop stress is present and
  contributes to the von Mises combination, but is roughly a third as large.
- The radial stress at the free edge is zero, and the hoop stress is small but not zero.
  The free-edge condition M_r = 0 is satisfied exactly; the hoop moment is not zero there,
  which is correct — a free edge releases radial moment, not hoop moment.
- The deflection is **2.1012 µm**, occurring at the free rim, and is monotonically
  increasing from zero at the clamp — again exactly as expected for a clamped plate.

### 7.1 The small-deflection assumption

The maximum deflection is 2.10 µm against a 6 mm thickness. That is

    2.10 µm / 6000 µm = 0.00035 = 0.035% of the thickness

so the geometry is changed by essentially nothing, and the linear small-deflection theory
is valid with an enormous margin. The project notes quote "0.09 percent", which is
**inconsistent** with the solver output and the deflection ratio reported in the
optimisation table (both of which give 0.088% for the *6 mm plate at the 3152 N ultimate
load*, and 0.035% at 1261 N). The discrepancy is one of bookkeeping between load cases,
not of physics; see [Section 12](#12-discrepancies-and-issues-found-in-the-written-record).

### 7.2 The stress distribution

![Through-width bending stress and transverse deflection at 1261 N](stress_profile.png)

The upper panel shows radial, hoop and von Mises stress across the width, with the yield
stress and the load-band edge marked. Note how far below the 250 MPa yield line the curves
sit: the plate is operating at about 2.5% of yield. The lower panel shows the deflection,
rising smoothly from zero at the clamp to the rim.

### 7.3 The beam comparison

| Quantity | Beam model | Plate model | Ratio |
| :---- | :---- | :---- | :---- |
| Peak bending stress | 161.67 MPa | 6.19 MPa | **26.1×** |
| Safety factor | 1.55 | 40.39 | 26.1× |

The solo purpose of `beam.py` is to make this factor visible and permanent.

---

## 8. Validation of the solver

**ANSYS results are not yet available.** No `ANSYS_RESULTS.txt` exists in the workspace —
this was verified directly. Therefore the solver has been validated only against
independent analytical and physical checks. All **36 tests** pass.

The validation instruments, and what each one actually proves:

| Check | Agreement | What it proves |
| :---- | :---- | :---- |
| Exact symbolic annulus solution, uniform load | 1e-9 relative | The basis functions, derivatives and assembly reproduce a known textbook closed form |
| Free edge M_r = 0 and Q = 0 | exact to roundoff | The free-edge conditions are actually satisfied by the solution |
| Clamped edge w = 0 and w' = 0 | exact to roundoff | The clamp conditions are actually satisfied |
| Clamp reaction equals applied load | 1e-5 relative | **Global equilibrium** — the total force is conserved through the model |
| Ring-of-force statics for clamp moment | 2.15% | An independent hand calculation of the clamp moment agrees to ~2% |
| Linear scaling of stress with load | 1e-9 relative | Linearity holds, so superposition and the √N scaling law are legitimate |
| Linear scaling of deflection with load | 1e-9 relative | Same |
| Cubic scaling of deflection with thickness | 1e-6 relative | The t³ stiffness law is correctly implemented, confirming D is right |
| Region splitting does not change the answer | 1e-9 relative | The continuity conditions are correct — splitting a uniform load into two regions changes nothing |
| Zero load gives zero result | exact | No spurious forcing or offset |

**The two most valuable checks.** Most of the above test the solver against itself. Two do
not:

1. **The exact symbolic annulus solution.** A separate, independent closed-form solution
   exists for a uniform load over a full annulus. The solver returns M_r = −25.0060 N·m/m
   and M_t = −8.2520 N·m/m at the clamp, matching the symbolic values to 1 part in 10⁹.
2. **The ring-of-force statics estimate.** This is a deliberately crude hand calculation
   that treats the distributed thrust as a single line load at the epoxy band's centroid
   radius (42.5 mm) and computes the clamp moment by simple lever arm:

       M = P (r_c − b) / (2πb) = 1261 × (0.0425 − 0.035) / (2π × 0.035) = 43.0059 N·m/m

   The solver gives −42.0799 N·m/m. That is a **2.15%** difference **in magnitude**. Note
   that the two carry opposite signs, which is expected: the hand estimate is a
   magnitude-of-moment calculation with no sign convention attached, while the solver
   reports the signed moment consistent with its curvature relation. Only the magnitudes
   are being compared here. Given that the hand
   estimate ignores the actual distribution of load across the band and the finite bending
   stiffness of the annulus, agreement to ~2% is exactly what one would hope for. This is
   the strongest *physical* — as opposed to mathematical — corroboration available.

**What this validation does and does not establish.** It establishes that the code
correctly and consistently solves the equations it was given, that those equations satisfy
global equilibrium, and that the implementation matches known exact solutions. It does
**not** establish that Kirchhoff thin-plate theory, with a rigid clamp and a uniform epoxy
pressure, faithfully represents this real part. Only ANSYS (or a test) can establish that.

---

## 9. Optimisation findings

### 9.1 The target

The requirement: safety factor 2.5 against yield at the 2.5× design load of 3152 N.
Therefore the allowable peak stress is

    250 MPa / 2.5 = 100 MPa

### 9.2 The closed-form result

Using the scaling law t ∝ √N:

    t_required = 6 · √(15.4741 / 100) = 2.3602 mm

verified to give a safety factor of 2.500000000000008 against the solver — the closed form
and the solver agree to about 1.6 parts in 10⁶ in thickness, and the safety factor lands on
2.500 to better than one part in 10¹⁴.

### 9.3 Comparison of designs

| Quantity | 6 mm as built | 2.360 mm optimised |
| :---- | :---- | :---- |
| Thickness | 6.0000 mm | 2.3602 mm |
| Peak von Mises at 3152 N | 15.4741 MPa | 100.0020 MPa |
| Safety factor at 3152 N | 16.16 | **2.50** |
| Peak von Mises at 1261 N | 6.1896 MPa | 40.00 MPa |
| Safety factor at 1261 N | 40.39 | 6.25 |
| Maximum deflection at 3152 N | 5.2531 µm | 86.3018 µm |
| Annulus mass | 111.1521 g | 43.7235 g |

    mass saved  = 67.4286 g
    mass saving = 60.66 %

![Von Mises stress at ultimate load and the mass comparison](optimisation.png)

### 9.4 Is the thinner plate still a valid thin plate?

Making the plate thinner makes the thin-plate assumptions harder to satisfy, so the
dimensionless ratios were re-checked. This is a genuinely important step and was not
skipped.

| Ratio | 6 mm | 2.360 mm | Requirement | Verdict |
| :---- | :---- | :---- | :---- | :---- |
| Deflection / thickness | 0.00088 | 0.03657 | ≪ 1 | Satisfied with room |
| Overhang / thickness | 3.33 | 8.47 | > ~5 | **Improved**; 6 mm actually fails this |
| Outer radius / thickness | 9.17 | 23.30 | > ~5 | Improved |
| Matrix condition number | 2.36e6 | 1.43e5 | not too large | Improved |

Two observations:

- **Both directional ratios improve** at the thinner section, because a thinner plate is
  *relatively* wider and thinner plates are exactly where thin-plate theory is most
  applicable. It is a pleasing result that the structurally efficient design is also the
  one the theory describes best.
- **The condition number improves** by more than an order of magnitude, so there is no
  numerical penalty. Both values are comfortably small; there is no conditioning concern at
  either thickness.
- The 6 mm plate's overhang/thickness of 3.33 is marginally below the "greater than about 5"
  guideline, meaning the *original* part is arguably the one that better needs a thick-plate
  or shear-corrected treatment. This is worth noting when comparing against ANSYS.

The deflection ratio rises from 0.09% to 3.7% of thickness but remains far inside the small
deflection regime. **The design is valid under the current model.**

### 9.5 The dominant risk in accepting 2.36 mm

This deserves its own emphasis. The original 6 mm plate had a safety factor of 16 at the
design load. Every modelling error therefore had a large cushion: even if the true stress
were 50% higher than predicted, the plate would still be safe by a factor of 10.

At 2.360 mm and a safety factor of 2.5, that cushion is gone. The known omissions —
transverse shear deformation, non-uniform epoxy bond stress, bracket compliance, the cable
hole stress concentration — would each **increase** the true stress. Any one of them
producing more than a ~2.5× error would make the design marginal, and two acting together
easily could. The 2.36 mm figure should therefore be regarded as the output of a model that
has passed all available self-consistency checks but whose physics remains unconfirmed.

### 9.6 The stress distribution and the case for a taper

At 2.360 mm and 3152 N, the stress across the width is strongly non-uniform:

| Radius, mm | von Mises, MPa | Fraction of peak |
| :---- | :---- | :---- |
| 35.0 | 100.00 | 100.0% |
| 38.0 | 59.06 | 59.1% |
| 40.0 | 39.04 | 39.0% |
| 42.0 | 24.16 | 24.2% |
| 45.0 | 11.34 | 11.3% |
| 48.0 | 8.85 | 8.8% |
| 50.0 | 8.60 | 8.6% |
| 55.0 | 7.41 | 7.4% |

Only material near the clamp is working hard. Beyond radius 45 mm the stress is already
below 12% of peak. **A uniform plate is therefore structurally wasteful.**

**The estimate.** Scaling local thickness as 1/√σ gives a profile falling from 2.360 mm at
the clamp to 0.642 mm at the rim. Applied to the 35–55 mm region, this would reduce that
region's mass from 36.0 g to about 15.8 g — a further 56% of that region. Over the whole
plate, including the 29–35 mm region under the bracket (7.7 g, which must stay at full
thickness to carry the clamp bearing), the tapered plate would weigh about **23.5 g against
43.7 g** for the uniform plate — a further **46%**.

**Why this estimate cannot be trusted yet.** The estimate proceeds by asking, for each
radius, "how thin could this local section be if it alone had to reach 100 MPa?" That is a
*pointwise* argument, and the solver is not a pointwise model. Two consequences:

1. **The solver assumes a constant flexural rigidity.** A variable thickness means D varies
   with r, which changes the governing equation and adds continuity conditions. The
   closed-form basis used here does not apply to a taper at all, so the profile cannot be
   verified — or iterated — without extending the solver. This is real work, not a
   parameter change.
2. **The pointwise profile is not self-consistent.** Making a section thinner raises the
   stress in it (σ ∝ 1/t²), which in turn demands a thicker section. The 1/√σ profile is
   therefore the *first* step of an iterative process, not its answer. A proper solution
   would iterate the thickness field against the recomputed stress until both agree — and
   each iteration changes the plate's stiffness distribution, which changes the stress
   everywhere, including at the clamp. The true optimum is therefore thicker than this
   estimate at the rim.
   For reference, the rim as estimated (0.6424 mm) would see a local deflection of 86.30 µm,
   a ratio of 0.134 to its own thickness — inside small-deflection theory in isolation, so
   large deflection is not the immediate obstacle it might appear. The obstacle is the
   iteration above.

The estimate should therefore be read as an **upper bound on the achievable saving**, not a
prediction. A real tapered design would land **somewhere between 23.5 g and 43.7 g**, and
needs the solver extended to variable thickness before any figure can be trusted.

---

## 10. Limitations

These are the known weaknesses of the analysis. They are grouped by severity, and each is
stated with its likely direction of error — because "this is neglected" is far less useful
than "this is neglected, and it would make the answer worse".

### 10.1 Blocking limitation

- **The physics is unvalidated.** No ANSYS results exist. The solver is validated against
  analytical and equilibrium checks only. Everything in this report rests on Kirchhoff
  thin-plate theory with a rigid clamp and uniform epoxy pressure being a faithful model.
  **All quantitative results should be treated as provisional until ANSYS confirms them.**

### 10.2 Known omissions that would increase the true stress

- **The 8 mm cable hole is not modelled.** *Effect: increases local stress.* The hole is a
  local feature that concentrates stress, so it does not change the global load path — the
  correct treatment is to apply a stress concentration factor (typically 2–3 for a plain
  circular hole in a bending field) to the local smooth-theory stress, rather than model the
  hole in the solver. This decision is still open. This is arguably the **largest unquantified
  item**, and at a target safety factor of 2.5 a factor of 3 would be fatal to the 2.36 mm
  design.
- **The epoxy bond stress is assumed uniform over the band.** *Effect: could increase or
  decrease peak stress, but peaks are more likely than troughs.* Real adhesive joints rarely
  distribute load evenly; they tend to concentrate near the edges of the bond, particularly
  the inner edge — which is exactly where the plate's peak stress already is. This could be
  a significant under-prediction and is not quantifiable without bond-line analysis.
- **Transverse shear deformation is neglected.** *Effect: increases deflection, modestly
  increases stress.* Kirchhoff theory assumes normals stay normal; a thick plate shears,
  which adds deflection beyond the bending prediction. With the 6 mm plate's
  overhang/thickness of 3.33 (below the ~5 guideline) this is a second-order but real effect.
  At 2.360 mm the ratio improves to 8.47 and the effect shrinks. It should nonetheless be
  confirmed — a Mindlin/Reissner thick-plate or a shear correction would settle it.

### 10.3 Assumptions that are probably conservative or neutral

- **The bracket is treated as perfectly rigid.** *Effect: reduces the true clamp moment.* A
  real bracket has finite stiffness and will rotate slightly under load, which relieves the
  clamped edge and lowers the peak moment. This assumption is therefore **conservative**:
  the real part is probably slightly less stressed than predicted.
- **Body forces and plate self-weight are neglected.** *Effect: negligible.* The plate weighs
  111 g against a 1261 N thrust — a ratio of under 0.1%. Agreed with the supervisor. During
  the 2.7 s burn the axial acceleration would add a small inertial load, but it is a
  small correction to a stress that is already only 2.5% of yield.

### 10.4 Optimisation-specific limitations

- **Only thickness is optimised.** The radii, load band and bracket diameter are all
  unchanged. The optimisation is one-dimensional.
- **Deflection is not constrained.** At 86 µm the optimised plate deflects 3.7% of its
  thickness. Structurally harmless, but it may matter functionally — a thrust plate that
  moves could affect motor alignment or the epoxy bond line. This has not been checked
  against any functional requirement.
- **Buckling is not checked.** A thin section under compressive thrust could in principle
  buckle. With a 20 mm span and 2.36 mm thickness this is not expected to govern, but it has
  **not been verified** — and at the tapered estimate of 0.64 mm it would become a serious
  concern.
- **Manufacturability is unconfirmed.** Whether 2.36 mm is practical to make, and whether
  the bracket and epoxy interfaces require a minimum bearing thickness, are both open
  questions. A taper to 0.64 mm would raise these questions sharply.

### 10.5 Method limitations

- **The model is strictly axisymmetric.** It cannot represent any non-axisymmetric effect:
  a locally stiffened bracket, an asymmetric bond line, bolt preload, or the local effect of
  the cable hole. All of these are inherently three-dimensional.
- **Boundary conditions are idealised.** A perfect clamp (zero rotation) and a perfect free
  edge are mathematical abstractions. Real structures lie between them. The clamp idealisation
  is the more consequential of the two.
- **Results are only as good as the input data.** Yield stress, modulus, band radii and
  bracket diameter are all taken as given. No tolerance or variability analysis has been
  performed.

---

## 11. What is still to be done

Ranked by importance:

1. **Validate against ANSYS.** This is the gating item. Until it is done, nothing else
   matters much. When `ANSYS_RESULTS.txt` appears, the comparison should be against the
   6 mm plate at both 1261 N and 3152 N, checking peak von Mises stress, its *location*, and
   peak deflection. The location check is the most diagnostic: if ANSYS also peaks at the
   clamp radius, the boundary condition modelling is confirmed; if it peaks elsewhere, the
   clamp or load idealisation is wrong.
2. **Decide the cable hole stress concentration factor**, and apply it. This may well
   change the optimal thickness materially.
3. **Extend the solver to variable thickness** if the taper is to be pursued, and accept
   that a large-deflection analysis will be needed at the thin end.
4. **Confirm manufacturability** of the chosen thickness, and check the bracket and epoxy
   interfaces for minimum bearing thickness requirements.
5. **Check buckling** of the thin section under thrust, at least by a standard closed-form
   plate buckling estimate.
6. **Consider a shear-deformation correction** (Mindlin/Reissner) as a cross-check,
   especially for the original 6 mm plate whose overhang/thickness of 3.33 is on the low
   side.
7. **Reconcile the 441.3 g weight anomaly** — see below.

---

## 12. Discrepancies and issues found in the written record

These were discovered while preparing this report and are recorded so they are not lost.
None affects the validity of the solver; all affect the documentation.

**1. The 6 mm plate's mass discrepancy (already identified, unresolved).** The 441.3 g
weight originally tabulated for the part is not consistent with a flat 26 mm wide, 6 mm
thick aluminium annulus, which weighs **111.152 g** — a factor of about 4 too heavy. The
user removed the figure from the tables, and it is correctly *not* used as the optimisation
baseline. However it has never been explained. Possible causes: the 441.3 g may refer to an
assembly (bracket plus plate), or to a different material, or may simply be an error. **This
is worth resolving**, because if the real part is heavier than 111 g for a reason not
captured by the annulus geometry, that reason may also affect the structural model.

**2. Test count is wrong in `CONTEXT.md`.** `CONTEXT.md` line 143 states "All 27 tests in
`test_plate.py` pass." The file actually contains **36** tests, all of which pass. The
figure appears to be a stale count from before tests were added. `README.md` correctly
avoids quoting a number.

**3. Deflection-as-fraction-of-thickness is quoted inconsistently.** `CONTEXT.md` says the
6 mm plate's "peak deflection is 0.035 percent of the plate thickness", while the
optimisation table in the same file gives 0.00088 for the 6 mm plate. Both are correct but
for **different load cases**: 0.035% is at 1261 N, 0.088% is at 3152 N. The README repeats
the 0.09% figure. The text does not make the load case explicit, so the two figures appear
to contradict each other. It should state the load case each time.

**4. `README.md` equilibrium tolerance wording.** The README's validation table lists the
symbolic annulus check as "Mr, Mt at clamp" agreeing to 1e-9. This is correct and was
re-verified (M_r = −25.0060247318874, M_t = −8.251988161522812 N·m/m). No issue — recorded
only to confirm it was checked.

**5. `CONTEXT.md` overstates the deflection at the tapered rim.** `CONTEXT.md` line 273
states: "at 0.64 mm the local deflection reaches about 6.8 times the local thickness, which
breaks the small deflection assumption." The solver gives a rim deflection of **86.30 µm**
against a rim thickness of 0.6424 mm — a ratio of **0.134**, i.e. about 7.4%, not 6.8. The
figure appears to be the *deflection-to-thickness ratio of the 6 mm plate's overhang*
(3.33) or a similar ratio, misapplied. The conclusion that the taper needs a variable-
thickness and possibly large-deflection treatment is still correct, but for the reason given
in [Section 9.6](#96-the-stress-distribution-and-the-case-for-a-taper) — the solver cannot
iterate a thickness field — rather than because the small-deflection limit is breached. This
report does not repeat the 6.8 figure.

**6. Minor: `CONTEXT.md` header structure.** The file uses both `#` and `##` top-level
headings inconsistently (e.g. "# Technical definition", "# Optimisation findings") where
the rest of the document is hierarchically nested. Cosmetic only.

**7. Note on `beam.py` being retained.** `beam.py` contributes nothing to the results but is
imported and exercised by `test_plate.py` (two tests depend on it). It is correctly labelled
as a comparison baseline in both the README and `CONTEXT.md`. This is good practice — it
preserves the record of a wrong turn. However, its presence means the test count includes
tests of a model that is deliberately not used, which could confuse a future reader. A brief
note in the test file would help.

---

## 13. Assessment and recommendation

**What has been achieved, stated fairly.** A complete, self-consistent, closed-form solver
for the stress and deflection of this annular plate has been built and extensively checked
against independent analytical and equilibrium criteria. It has been used to identify that
the existing part is dramatically over-designed, and to compute the thickness that would
meet the stated safety requirement — a 60.7% mass saving in closed form, with no search.
The identification of the taper opportunity, though not yet quantifiable, is a real and
potentially larger finding. The work is methodologically sound: the incorrect beam model was
identified, understood, corrected, and retained as a documented comparison rather than
quietly deleted.

**What has not been achieved.** The physics has not been validated. The CV deliverable
explicitly requires agreement "within [X]% of ANSYS", and that X cannot be filled in yet.
Until it is, the solver is a well-tested hypothesis, not a validated result.

**The recommendation.** Treat **2.3602 mm as the model's prediction, pending validation**.
Do not present it as an accepted design. The most valuable next step by a wide margin is the
ANSYS comparison, and the second most valuable is deciding the cable-hole stress
concentration treatment — because at a safety factor of 2.5 there is very little margin, and
a stress concentration factor of 2–3 applied on top would invalidate the design outright.

**A caution about the headline saving.** The 60.7% figure is correct relative to a 6 mm flat
plate of the stated geometry. It is a large relative saving on a very small component: in
absolute terms, 67 g on a roughly 11 kg rocket, or about 0.6% of the vehicle mass. That is
worth stating plainly so the result is not oversold. The engineering value of the work lies
more in the validated method and the understanding of the load path than in the grams saved
— though if the taper is successfully pursued, the saving roughly doubles.

---

## Appendix A — Quick-reference data

**Geometry**

| Quantity | Symbol | Value |
| :---- | :---- | :---- |
| Outer radius | r_o | 55.0 mm |
| Inner radius | r_i | 29.0 mm |
| Bracket (clamp) radius | b | 35.0 mm |
| Unsupported overhang | r_o − b | 20.0 mm |
| Ring width | r_o − r_i | 26.0 mm |
| Load band | — | 35.0 – 50.0 mm |
| Thickness (baseline) | t | 6.0 mm |
| Thickness (optimised) | t* | 2.3602 mm |

**Material — Al 6082-T6**

| Quantity | Symbol | Value |
| :---- | :---- | :---- |
| Young's modulus | E | 70 GPa |
| Poisson's ratio | ν | 0.33 |
| Density | ρ | 2700 kg/m³ |
| Yield stress | σ_y | 250 MPa |
| Ultimate stress | σ_u | 300 MPa |

**Loads**

| Quantity | Value |
| :---- | :---- |
| Peak thrust | 1261 N |
| Design (ultimate) thrust | 3152 N |
| Epoxy band area | 4005.53 mm² |
| Pressure at 1261 N | 314.81 kPa |
| Pressure at 3152 N | 786.91 kPa |

**Key results**

| Quantity | 1261 N | 3152 N |
| :---- | :---- | :---- |
| Peak von Mises, 6 mm | 6.1896 MPa | 15.4741 MPa |
| Peak von Mises, 2.3602 mm | 40.00 MPa | 100.00 MPa |
| Deflection, 6 mm | 2.1012 µm | 5.2531 µm |
| Deflection, 2.3602 mm | 34.5196 µm | 86.2991 µm |
| Safety factor (6 mm) | 40.39 | 16.16 |
| Safety factor (2.3602 mm) | 6.25 | 2.50 |

**Masses**

| Configuration | Mass |
| :---- | :---- |
| 6.0000 mm annulus | 111.1521 g |
| 2.3602 mm annulus | 43.7235 g |
| Saving | 67.4286 g (60.66%) |
| Tapered estimate (indicative only) | ~23.5 g |

## Appendix B — Reproducing these results

    python3 run.py
    .venv/bin/python -m pytest test_plate.py -q

The virtual environment is created with `--system-site-packages` so that NumPy and SymPy
remain available:

    python3 -m venv --system-site-packages .venv
    .venv/bin/pip install pytest

`run.py` prints the full summary and writes `stress_profile.png` and `optimisation.png`.
Every number in this report was reproduced by running the code during its preparation. The
suite reports **36 passed**.
