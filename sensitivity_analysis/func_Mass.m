function [mass] = func_Mass(vals, xmesh)
    mass = trapz(xmesh, vals/max(vals));
end