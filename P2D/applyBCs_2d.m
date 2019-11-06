% ---------------------------------------------------------------
% retrieve boundary
% ---------------------------------------------------------------
nodes_left = bound.left;
nodes_right = bound.right;
% left boundary of cathode
node_caleft = bound.caleft;

if itime == 1

    % in this case, no interfacial flux
    % this is already ensured
    % K(dofsFlux, :) = 0;
    % K(dofsFlux, dofsFlux) = eye(length(dofsFlux), length(dofsFlux));
    % Res(dofsFlux) = -(0 - u(dofsFlux, itime));

    % then you need to specify a reference potential for electrolyte
    rw = dofArray(nodes_left, 1);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(phi0_el - u(rw, itime));

    % reference potential for anode
    rw = dofArray(nodes_left, 3);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(phi0_an - u(rw, itime));

    % reference potential for cathode
    rw = dofArray(nodes_right, 3);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(phi0_ca - u(rw, itime));

else % for other time steps

    % ---------------------------------------------------------------
    % boundary condition for anode and cathode
    % ---------------------------------------------------------------

    % current for electrode material
    % constant current flux at the right edge of electrode
    Fr = applyCurrent(ncoord, coords, ndofs, nfacenodes, dofArray, coefs, bound);

    Res = Res + Fr;

    % fixed potential for electrode at anode left
    rw = dofArray(nodes_left, 3);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(0 - u(rw, itime));

end
