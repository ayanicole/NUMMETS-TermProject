# NUMMETS-TermProject
This is our group's workspace for the MATLAB code, figures, and IEEE paper. Everything we do goes here so we always have one up-to-date version.

## Model Lead modules (grid, assembly, heat flow, interpolation)

Run `test_model` from the repo folder. The last line should read `ALL CHECKS PASSED` (28 checks, tested in GNU Octave 10.3).

| File | What it does |
|---|---|
| `build_grid.m` | Grid and indexing. Node (i,j) is unknown p = i + (j-1)*Nx. `reshape(T,Nx,Ny).'` gives the Ny x Nx map. |
| `assemble_system.m` | `[A,b,nodes]`: the linear 5-point system, with half/quarter-cell energy balances on edges and corners. Radiation is returned separately in `nodes.rho`. |
| `radiative_bc.m` | Newton-Raphson for one radiating node: `radiative_bc(Tguess,P,tolBC,Tnb,beta)`. The old 3-argument stub could not work, because the node equation also needs its neighbours and weight. |
| `heat_flow.m` | Edge heat flow by one-sided 3-point flux plus trapezoid or Simpson. `'radiation'` mode closes the energy balance. |
| `make_source.m` | Uniform, rectangular or Gaussian heat source with exact total power. Always build sources with this. |
| `material_k.m` | Newton divided-difference k(T) for copper and aluminum. |
| `interp_field.m` | Lagrange interpolation of the temperature field. |
| `test_model.m` | The 28 checks. |

**Full nonlinear system:** `F(T) = A*T - b + nodes.rho.*(T.^4 - Tamb^4) = 0`.

**For the Solver Lead:** for a radiating node p inside a sweep, use
`d=A(p,p); Tnb=(b(p)-A(p,:)*T+d*T(p))/d; T(p)=radiative_bc(T(p),P,tolBC,Tnb,nodes.rho(p)/d);`
See `test_model` check 5c. With every edge radiating or insulated and little top-face loss, Gauss-Seidel needs thousands of sweeps, so keep one Dirichlet edge or use h >= ~100 W/m^2-K.

**For the Integrator (`parameters.m`):** the new fields are optional and off by default. They are `P.bc.S/N/W/E` (type and value), `P.q`, `P.h` and `P.t`. The defaults reproduce the current setup.

**Open findings for the group**
1. The defaults give a flat map (spread under 0.1 K), so copper and aluminum look identical. Use a source with `h` and `t` set.
2. With 10 W, t = 2.5 mm, h = 1000 and a 4 mm hotspot, copper peaks at 324.9 K and aluminum at 326.3 K.
3. Edge radiation is small at realistic temperatures (0.3% of the power at h = 1000). Use a low-h case to show the nonlinear solver.
