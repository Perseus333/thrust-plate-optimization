# Rocket Thrust-plate Analytical Optimisation

Closed-form stress and deflection solver for the annular thrust plate of a high-power
rocket, used to size the plate for a required safety factor.

## The problem

A flat aluminium annulus (Al 6082-T6) is loaded on its inner annulus by motor thrust
transmitted through an epoxy band, and is clamped by epoxy over its outer band. The
supported inner and outer edges are free.

| Property | Value |
| :---- | :---- |
| Plate outer diameter | 110 mm |
| Bore diameter | 58 mm |
| Thickness as built | 6 mm |
| Load annulus | radius 29 - 36 mm |
| Free plate | radius 36 - 43 mm |
| Clamped band | radius 43 - 55 mm |
| Clamp radius used | 49 mm, middle of the clamped band |
| Load area | 1429.4 mm^2 |
| Peak thrust | 1261 N |
| Design thrust, 2.5 safety | 3152 N |
| Required safety factor | 2.5 against yield (250 MPa) |
| Allowable peak stress | 100 MPa |

Load acts normal to the top face only. There is no load or fixing on the cylindrical
edges. The grey region between the load and the clamp is plain full-thickness plate.

## Method

Axisymmetric Kirchhoff thin-plate bending. Radial symmetry reduces the problem to one
coordinate, r. The deflection satisfies

    D * laplacian(laplacian(w)) = q(r),      D = E t^3 / (12 (1 - nu^2))

with homogeneous solution

    w(r) = C1 + C2 r^2 + C3 log(r) + C4 r^2 log(r)

and particular integral for uniform band pressure

    w_p(r) = q r^4 / (64 D)

Moments, shear and stresses follow in closed form:

    M_r = -D ( w'' + nu w' / r )
    M_t = -D ( w' / r + nu w'' )
    Q   = -D d/dr ( w'' + w' / r )
    sigma_r = 6 M_r / t^2
    sigma_t = 6 M_t / t^2

The domain is split at each load discontinuity, and regions are joined by continuity of
w, w', M_r and Q. The four remaining boundary conditions are a free edge at the bore and
a clamped edge at the clamp radius:

    r = 29 mm (bore)   :  M_r = 0  and  Q = 0
    r = 49 mm (clamp)  :  w = 0    and  w' = 0

This gives a small dense linear system solved directly. No numerical differentiation is
used: for the homogeneous basis, d/dr ( w'' + w'/r ) = 4 C4 / r.

The green band is a finite-area restrained annulus (43 - 55 mm), not a point clamp. It is
idealised as a single rigid clamp at its middle radius, r = 49 mm. Applying the clamp at
the two band edges was checked and shown to bracket the answer, so the private case in
which the plate is prevented from rotating at all under the band is central between them
and a single result is reported.

## Results

At the peak thrust of 1261 N, 6 mm plate:

| Quantity | Value |
| :---- | :---- |
| Peak von Mises stress | 9.78 MPa |
| Peak location | r = 49 mm, the clamp |
| Maximum deflection | 6.35 um |
| Safety factor against yield | 25.6 |

Optimised for a safety factor of 2.5 at the 3152 N design load:

| Quantity | Value |
| :---- | :---- |
| Required thickness | 2.967 mm |
| Peak von Mises at 3152 N | 100.00 MPa |
| Achieved safety factor | 2.50 |
| Deflection at 3152 N | 131.4 um |
| Plate mass | 54.96 g |
| Mass saving against 111.15 g | 50.6 % |

Because peak bending stress scales as 1 / t^2 at fixed load, the required thickness
follows in closed form without any search:

    t_required = t_reference * sqrt( sigma_reference / (yield / N) )

## Secondary checks

Neither transverse shear nor buckling governs.

Transverse shear, at the 3152 N load: the free span over thickness is 2.2 at 6 mm and 4.4
at 2.97 mm, low enough to ask whether thick-plate shear matters. It does not. The peak
transverse shear force is 13.9 N/mm at the load edge, a mean shear stress of 7.05 MPa at
the optimised thickness, which adds 0.7 % to the combined von Mises stress and 1.1 % to
the deflection. The shear capacity at the optimised thickness is 63 kN, twenty times the
load.

Buckling, at the 3152 N load: the annulus critical pressure is 29.6 MPa against an applied
2.21 MPa, a factor of 13.4, and the 7 mm free length is far too short to act as a column.
The narrowest margin is the hoop compression under the load band, where the classical ring
buckling force is 37.0 kN against the 15.4 kN carried, a factor of 2.4. The thin inner
section is therefore not at risk of elastic instability, but its buckling margin is of the
same order as its bending margin and should be revisited if the thickness is reduced
further.


