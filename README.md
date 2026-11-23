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
a clamped edge at the outer band:

    r = 29 mm (bore)  :  M_r = 0  and  Q = 0
    r = 43 mm (clamp) :  w = 0    and  w' = 0

This gives a small dense linear system solved directly. No numerical differentiation is
used: for the homogeneous basis, d/dr ( w'' + w'/r ) = 4 C4 / r.

The green band is a finite-area restrained annulus (43 - 55 mm), not a point clamp. Both
limiting idealisations are reported: a clamp at r = 43 mm (band start) and a clamp at
r = 55 mm (outer rim). The truth lies between them.

## Results

At the peak thrust of 1261 N, 6 mm plate:

| Quantity | Clamp at 43 mm | Clamp at 55 mm |
| :---- | :---- | :---- |
| Peak von Mises stress | 7.17 MPa | 11.65 MPa |
| Peak location | r = 43 mm | r = 55 mm |
| Maximum deflection | 2.10 um | 13.41 um |
| Safety factor against yield | 34.9 | 21.5 |

Optimised for a safety factor of 2.5 at the 3152 N design load:

| Quantity | Clamp at 43 mm | Clamp at 55 mm |
| :---- | :---- | :---- |
| Required thickness | 2.540 mm | 3.239 mm |
| Peak von Mises at 3152 N | 100.00 MPa | 100.00 MPa |
| Achieved safety factor | 2.50 | 2.50 |
| Deflection at 3152 N | 69.3 um | 213.2 um |
| Plate mass | 47.05 g | 60.00 g |
| Mass saving against 111.15 g | 57.7 % | 46.0 % |

Because peak bending stress scales as 1 / t^2 at fixed load, the required thickness
follows in closed form without any search:

    t_required = t_reference * sqrt( sigma_reference / (yield / N) )

