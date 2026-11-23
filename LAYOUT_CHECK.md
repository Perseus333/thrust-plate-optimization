# Corrected Layout Analysis — model vs. the real part

**Verdict: the original solver modelled the part inside-out.** It was fixed to a radius
*inboard* of the load and left the outer rim free. The real part is loaded at the **inner**
radius and fixed at the **outer** rim. The corrected model has been built, validated, and
run. The optimisation result changes.

---

## 0. Correction of an earlier error in this document

An earlier version measured the top-view image in pixels and inferred green = 48.6–55 mm,
red = 28.4–35.9 mm. **The image is not to scale, as confirmed by the user. Those numbers are
withdrawn.** Scaling a hand-drawn schematic and building a conclusion on it was a mistake.
The real radii, given by the user, are used everywhere below.

---

## 1. The real geometry (user-supplied, authoritative)

| Feature | Radial extent | Condition |
| :---- | :---- | :---- |
| Central hole | 0 → 29 mm | Empty |
| **Red — load annulus** | **29 → 36 mm** | Thrust on the **top face only**, perpendicular |
| **Grey — plain plate** | **36 → 43 mm** | Free, full thickness |
| **Green — clamped annulus** | **43 → 55 mm** | Fixed by epoxy on **both faces** |
| Outer rim | 55 mm | No free overhang — plate ends where the clamp ends |

Additional constraints stated by the user:

- The grey span is plain full-thickness plate.
- The green clamp is an annulus with **no free overhang**: the plate's outer radius is also
  55 mm.
- The hole boundary coincides with the red band's inner edge (29 mm), just as the green
  band's outer boundary coincides with the plate's outer edge (55 mm).
- **There are no loads or fixings on the cylindrical (vertical) edges** — only on the flat
  top and bottom faces.

Radial order, centre outward:

    HOLE  →  LOAD  →  FREE  →  CLAMPED  →  outside

Load annulus area: π(36² − 29²) = **1429.42 mm²** → pressure at 1261 N = **882.17 kPa**.

---

## 2. What the original solver assumed

    HOLE → CLAMPED(35 mm) → LOAD BAND(35–50 mm) → FREE EDGE(55 mm)

The two orders are reversed. In the original model the load sat *outboard* of the support
so the plate cantilevered outward; in reality the load is *inboard* of the support, so the
plate is a centrally loaded disc held at its rim. The peak moment consequently belongs at
the **outer** edge, not the inner radius.

### The solver could not express the real layout

`plate.assemble` applies the clamp to region 0 (the innermost) and the free condition to the
last region, so `clamp_radius` was structurally forced to be the inner boundary. Two
consequences:

- A load band lying outside the swept domain was **silently dropped** by `collect_cuts`
  (`plate.py` lines 85–92), producing a uniformly zero solution that looked valid. That is a
  robustness defect worth fixing independently of this layout question.
- Passing an outer clamp radius misapplied the boundary conditions. Verified: with
  clamp = 43 mm the "solution" had w(43) = −5.8e−08 and M_r(29) = −24.6 instead of 0.

The **mathematics** was never the obstacle — the basis functions, derivative relations and
assembly satisfy an outer-clamp problem to machine precision. Only the boundary-condition
placement and region-building logic were rigid.

---

## 3. The corrected model

New files (the originals are left untouched for comparison):

| File | Purpose |
| :---- | :---- |
| `plate_corrected.py` | Plate solver with the clamp placeable at **either** radial end |
| `plate_layout.py` | The real geometry and load case |
| `test_plate_layout.py` | 15 tests for the corrected layout — all passing |

The key change: boundary conditions are now attached to **named edges**
(`plate_corrected.Edge(radius, condition)`) instead of being inferred from position. A load
band outside the solved domain now raises `ValueError` instead of returning zeros.

### Validation of the corrected model

All 15 tests pass, including:

| Check | Result |
| :---- | :---- |
| Global equilibrium, clamp at 43 mm | reaction / applied = **−1.000000** |
| Global equilibrium, clamp at 55 mm | reaction / applied = **−1.000000** |
| Free edge M_r = 0 at the bore | 2.6e−13 (roundoff) |
| Clamped edge w = 0 at the rim | 0.0 exactly |
| Ring-of-force statics, clamp 43 mm | **0.6%** |
| Ring-of-force statics, clamp 55 mm | **3.5%** |
| Linear scaling of stress with load | 1e−9 relative |
| Cubic scaling of deflection with thickness | 1e−6 relative |
| Region splitting does not change the answer | 1e−9 relative |
| Out-of-domain load band raises | passes |

The independent ring-of-force check treats the thrust as a line load at the band centroid
(32.5 mm) giving a lever arm to the clamp — a completely separate calculation from the
solver, and it agrees to well under 4%.

---

## 4. An open modelling choice: which radius is clamped

The green band is 43 → 55 mm, fixed on both faces. It is a **finite-area restraint**, not a
point clamp, so two idealisations are defensible and both were run:

- **(A) Clamp at r = 43 mm.** The 43–55 mm annulus is treated as a rigid fixture; the
  bending plate is 29–43 mm. Conservative in the sense that it assumes the clamp gives no
  rotational compliance at all and ignores the bending of the 43–55 ring itself.
- **(B) Clamp at r = 55 mm.** The whole 29–55 mm is treated as plate, held only at the outer
  rim. These are bounds, and the truth lies between them — closer to (A), because the epoxy
  really does restrain the 43–55 mm band, but not perfectly rigidly.

**A full answer needs the clamp modelled as a distributed translational and rotational
spring over 43–55 mm**, which requires an extension beyond simple edge conditions. Until
then, (A) and (B) bracket the answer.

---

## 5. Corrected results

### 5.1 The 6 mm plate as built

| Quantity | Model (A) clamp 43 mm | Model (B) clamp 55 mm | *Old (wrong topology)* |
| :---- | :---- | :---- | :---- |
| Peak von Mises at 1261 N | **7.17 MPa** | **11.65 MPa** | *6.19 MPa* |
| Peak location | r = 43 mm | r = 55 mm | *r = 35 mm* |
| Peak radial stress | −8.12 MPa | −13.21 MPa | *−7.01 MPa* |
| Peak hoop stress | −2.68 MPa | +5.76 MPa | *−2.31 MPa* |
| Max deflection at 1261 N | 2.103 µm | 13.414 µm | *2.101 µm* |
| Safety factor at 1261 N | **34.9** | **21.5** | *40.4* |
| Peak von Mises at 3152 N | 17.92 MPa | 29.14 MPa | *15.47 MPa* |
| Safety factor at 3152 N | 13.95 | 8.58 | *16.2* |

Note that the peak stress is **higher** than the old model predicted in both cases (by 16%
and 88% respectively), and the deflections differ much more. The conclusion that the part is
heavily over-designed is unchanged — the safety factor is 21–35 rather than 40, still far
above any requirement.

### 5.2 The clamp location is not a detail

Models (A) and (B) differ by a factor of **1.63** in peak stress and **6.4** in deflection.
This confirms the point raised before the geometry was known: the answer is highly sensitive
to where the restraint is, and the restraint is the one thing the drawings describe
qualitatively but not precisely.

### 5.3 Optimised thickness for SF 2.5 at 3152 N

Since peak stress still scales as 1/t², the closed form
`t = 6·√(σ_ref / 100)` applies unchanged.

| | Model (A) clamp 43 mm | Model (B) clamp 55 mm | *Old (wrong topology)* |
| :---- | :---- | :---- | :---- |
| Required thickness | **2.540 mm** | **3.239 mm** | *2.360 mm* |
| Achieved peak at 3152 N | 100.00 MPa | 100.00 MPa | *100.00 MPa* |
| Achieved safety factor | 2.5000 | 2.5000 | *2.5000* |
| Deflection at 3152 N | 69.3 µm | 213.2 µm | *86.3 µm* |
| Plate mass | **47.06 g** | **60.00 g** | *43.72 g* |
| Mass saving vs 111.15 g | **57.7%** | **46.0%** | *60.7%* |

**The headline number moves from a 60.7% saving to between 46% and 58%**, with the more
physically representative idealisation (A) giving **57.7%** and a thickness of **2.54 mm**.

The 3.239 mm figure under model (B) sits closer to the original design intent and would be
the more conservative choice to carry forward if the clamp's true rotational stiffness is
unknown — but (B) is deliberately pessimistic, since it ignores the restraint over
43–55 mm entirely.

---

## 6. What this changes in the main report

`REPORT.md` documents the original model and its results. The following are now superseded:

| Claim in REPORT.md | Status |
| :---- | :---- |
| Clamp at r = 35 mm, free at r = 55 mm | **Wrong** — the clamp is the outer band, 43–55 mm |
| Load band 35–50 mm, area 4005.53 mm² | **Wrong** — the load annulus is 29–36 mm, area 1429.42 mm² |
| Pressure at 1261 N = 314.81 kPa | **Wrong** — it is 882.17 kPa |
| Peak von Mises 6.19 MPa at r = 35 mm | **Superseded** — 7.17 MPa (A) / 11.65 MPa (B), at the clamp |
| Deflection 2.10 µm | **Superseded** — 2.10 µm (A) / 13.41 µm (B) |
| Safety factor 40.4 | **Superseded** — 34.9 (A) / 21.5 (B) |
| Optimum thickness 2.360 mm | **Superseded** — 2.540 mm (A) / 3.239 mm (B) |
| Mass saving 60.7% | **Superseded** — 57.7% (A) / 46.0% (B) |

Still valid from the main report:

- The **method** (axisymmetric Kirchhoff plate, closed-form basis, region matching, no
  numerical differentiation) is correct and unchanged.
- The **beam is the wrong idealisation** conclusion — though the factor of 26 would change.
- The **scaling law** `t ∝ √N` and the closed-form optimiser are unchanged.
- The **limitations** (cable hole, epoxy non-uniformity, transverse shear, bracket
  compliance) all still apply, and the reduced safety margin argument still holds.
- The **validation approach** was sound; it simply validated the wrong geometry.

---

## 7. Remaining open items for this layout

1. **The 43–55 mm clamp: spring or rigid?** The largest remaining uncertainty. Modelling it
   as a distributed elastic restraint would replace the (A)/(B) bracket with one answer.
2. **Is the 8 mm cable hole within the clamped band or the free span?** This determines
   whether it sees significant stress at all — it may be far less critical than previously
   assumed, since the peak stress region is now the outer clamp.
3. **Transverse shear — now a first-order concern, not a second-order one.** The free span is
   only 36 → 43 mm = 7 mm. Span over thickness is therefore **1.17 at t = 6 mm** and
   **2.76 at t = 2.54 mm**, far below the ~5 guideline above which thin-plate theory is
   normally trusted. The full bending domain (29 → 43 mm = 14 mm) gives 2.33 and 5.51.
   These are *thick-plate* proportions: Kirchhoff theory neglects transverse shear, and at a
   span/thickness ratio near 1 that neglect is no longer small. **This is now the most
   likely source of a genuine modelling error**, and a Mindlin/Reissner or shear-corrected
   check should be treated as required rather than optional.
   For balance, the deflection-over-thickness ratios remain small — 0.00035 and 0.0273 under
   model (A), 0.0022 and 0.0658 under (B) — so the *small deflection* assumption is still
   comfortably satisfied. It is the *thin plate* assumption, not the small-deflection one,
   that is under strain.
4. **ANSYS validation**, still the gating item.
5. **The 441.3 g mass anomaly** remains unexplained.
