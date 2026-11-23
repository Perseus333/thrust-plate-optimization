# Rocket Thrust-plate Analytical Optimisation

> Stress and deflection solver for the annular thrust plate of a high-power rocket, used to size the thickness of the plate for the required safety factor.

## Context

A flat aluminium annulus (Al 6082-T6) is loaded on its inner top band by motor thrust transmitted through a bracket, and is clamped by epoxy over its outer band at the top and bottom. The inner and outer edges are free. Features on the plate such as holes are ignored for a simplified solution but taken into account within the 2.5x safety margin.

| Property | Value |
| :---- | :---- |
| Plate outer diameter | 110 mm |
| Bore diameter | 58 mm |
| Thickness as built | 6 mm |
| Load annulus | radius 29 - 36 mm |
| Clamped band | radius 43 - 55 mm |
| Clamp radius used | 49 mm |
| Peak thrust | 1261 N |
| Yield Stress | 250 MPa |
| Safety factor | 2.5 |
| Allowable peak stress | 100 MPa |

## Method



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

