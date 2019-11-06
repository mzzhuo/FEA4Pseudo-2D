F_rs = zeros(ndofs, 1);

node_right = bound.right;

rw = dofArray(node_right, 1);
% minus ensures: when pore wall flux is positive,
% species flux out of the particle
F_rs(rw) = -coefs.q_bc * pwflux;

if (itime == 1)
    F_rs = 0.0; % this means no flux out in the edge!
end

Res = Res + F_rs;
