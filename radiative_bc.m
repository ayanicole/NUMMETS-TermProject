function [Tnew, iters] = radiative_bc(Tguess, P, tolBC, Tnb, beta)
% RADIATIVE_BC  Solve one radiating boundary-node equation by Newton-Raphson.
%
%   F(T) = T - Tnb + beta*(T^4 - Tamb^4) = 0
%
%   where, for node p of the system built by assemble_system,
%       beta = nodes.rho(p) / A(p,p)
%       Tnb  = ( b(p) - A(p,:)*T + A(p,p)*T(p) ) / A(p,p)    % neighbours' contribution
%   (Tnb is what the node would be with no radiation, given its neighbours).
%
%   F is increasing and convex for T > 0 (F' = 1 + 4*beta*T^3, F'' = 12*beta*T^2),
%   so Newton-Raphson converges from any positive guess: monotonically from the
%   right of the root, and after one overshoot from the left.
%
% INPUT : Tguess (K), P (uses P.Tamb), tolBC (K, stop when |step| < tolBC),
%         Tnb (K), beta (1/K^3)
% OUTPUT: Tnew (K), iters (Newton steps used)
%
% NOTE: the original stub was radiative_bc(Tguess, P, tolBC). That cannot work
% because the node equation also needs its neighbours (Tnb) and its weight (beta).

    maxit = 50;
    T = Tguess;
    for iters = 1:maxit
        F  = T - Tnb + beta*(T^4 - P.Tamb^4);
        dF = 1 + 4*beta*T^3;
        dT = -F/dF;
        T  = T + dT;
        if abs(dT) < tolBC, break; end
    end
    if abs(dT) >= tolBC
        warning('radiative_bc:noconv', 'Newton did not reach tolBC in %d steps.', maxit);
    end
    Tnew = T;
end
