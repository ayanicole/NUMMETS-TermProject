# Model Lead hand-off

Tested in GNU Octave 10.3 (MATLAB itself was not available): run `test_model` from the repo folder. Expected last line: `ALL CHECKS PASSED`.

## Files
| File | Status | What it does |
|---|---|---|
| build_grid.m | replaces stub | Grid + indexing contract (node (i,j) -> p = i + (j-1)*Nx; `reshape(T,Nx,Ny).'` gives the Ny x Nx map) |
| assemble_system.m | replaces stub | `[A,b,nodes]`: linear 5-point system; radiation returned separately in `nodes.rho` |
| heat_flow.m | replaces stub | Edge heat flow, 'trap'/'simpson', optional edge and 'radiation' mode |
| radiative_bc.m | **changes signature** | `radiative_bc(Tguess,P,tolBC,Tnb,beta)`; the old 3-argument form lacks the neighbour data the node equation needs |
| make_source.m | new | Hotspot/uniform/Gaussian source with exact total power |
| material_k.m, interp_field.m | new | Newton divided-difference k(T); Lagrange probe of the field |
| test_model.m | new | 28 checks |

## For the Solver Lead
Full nonlinear system: `F(T) = A*T - b + nodes.rho.*(T.^4 - Tamb^4) = 0`.
For a radiating node p inside a sweep: `d=A(p,p); Tnb=(b(p)-A(p,:)*T+d*T(p))/d; T(p)=radiative_bc(T(p),P,tolBC,Tnb,nodes.rho(p)/d);` (non-radiating nodes: `T(p)=Tnb`). `test_model` check 5c runs exactly this loop.
Conditioning warning: with every edge radiating or insulated and little top-face loss, A is nearly singular and Gauss-Seidel needs thousands of sweeps. Keep a Dirichlet edge or h >= ~100 W/m^2-K.

## For the Integrator (parameters.m)
All new fields are optional and default to off: `P.bc.S/N/W/E` (type + value), `P.q` (use `make_source`), `P.h`, `P.t`. Defaults reproduce the current setup: south edge at Thot, other edges radiating.

## Findings the group must decide on
1. **Defaults give a flat map.** Hot south edge at 358 K, three radiating edges, no source: the whole plate is within 0.047 K of 358 K (copper) and about 0.08 K (aluminum). Copper vs aluminum would be indistinguishable.
2. **Meaningful scenario (10 W total, t = 2.5 mm, 20 mm plate, 41 x 41, h = 1000 W/m^2-K, all edges radiating):** uniform source spread 0.00 K; 4 mm square hotspot peak 324.9 K (Cu, spread 2.6 K) vs 326.3 K (Al, spread 4.2 K). Mean temperature (~323 K) is set by h, not by material.
3. **Radiation is small at realistic temperatures:** 0.3% of the power at h = 1000, 7% at h = 100 (mean 530 K), 75% at h = 10 (mean 928 K). Say so in the Limitations, and use a low-h case to show the nonlinear solver where it matters.
4. **Equal power:** a sharp-edged source sampled on nodes generated 12.7 W instead of 10 W. Always build sources with `make_source`.
5. **Energy balance:** on radiating edges use `heat_flow(...,'trap',edge,'radiation')`; sum of edge flows equals `nodes.Psrc` exactly (tested including corners).
