function k = material_k(name, T)
% MATERIAL_K  Thermal conductivity k(T) [W/(m*K)] of the spreader material,
% by Newton's divided-difference interpolation of tabulated values.
%
%   k = material_k('copper', 330)       % scalar or array T [K]
%
% TABLE VALUES: copper and aluminum (pure) read from Appendix A, Table A.1 of
% "Process Heat Transfer: Principles and Applications" (as retrieved online;
% the page layout was garbled). CONFIRM these against the textbook the group
% cites before the paper goes out, and update Td/kd below if they differ.
%
% Valid for 273 K <= T <= 800 K. Outside that range the polynomial would be
% extrapolating, so the function errors instead of guessing.

    switch lower(name)
        case 'copper'
            Td = [273 400 600 800];   kd = [401 392 383 371];
        case 'aluminum'
            Td = [273 400 600 800];   kd = [236 240 232 220];
        otherwise
            error('material_k:name', 'Unknown material "%s" (use ''copper'' or ''aluminum'').', name);
    end
    if any(T(:) < Td(1)) || any(T(:) > Td(end))
        error('material_k:range', 'T must be within %g..%g K.', Td(1), Td(end));
    end

    c = newton_coeffs(Td, kd);
    k = newton_eval(Td, c, T);
end

function c = newton_coeffs(x, y)
% Divided-difference coefficients c(1..n) of Newton's form.
    n = numel(x);
    c = y(:).';
    for j = 2:n
        for i = n:-1:j
            c(i) = (c(i) - c(i-1)) / (x(i) - x(i-j+1));
        end
    end
end

function v = newton_eval(x, c, xq)
% Nested (Horner-like) evaluation of Newton's polynomial.
    n = numel(x);
    v = c(n)*ones(size(xq));
    for i = n-1:-1:1
        v = v.*(xq - x(i)) + c(i);
    end
end
