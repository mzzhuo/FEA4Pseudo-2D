node_right = bound.right;
dofs_bound = dofArray(node_right, 1);
c_ss = u(dofs_bound,itime);

F_rs = zeros(ndofs, 1);
F_rs(dofs_bound) = -coefs.q_bc * 1;

deltaU = -K \ F_rs;

dcssdJ = deltaU(dofs_bound);
