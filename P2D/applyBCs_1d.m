
% ---------------------------------------------------------------
% retrieve boundary 
% ---------------------------------------------------------------
node_left = bound.left;
node_right = bound.right;
% left boundary of cathode
node_caleft = bound.caleft;

if itime == 1

    % in this case, no interfacial flux
    % this is already ensured
    % K(dofsFlux, :) = 0;
    % K(dofsFlux, dofsFlux) = eye(length(dofsFlux), length(dofsFlux));
    % Res(dofsFlux) = -(0 - u(dofsFlux, itime));

    % then you need to specify a reference potential for electrolyte
    rw = dofArray(node_left, 1);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(phi0_el - u(rw, itime));

    % reference potential for anode
    rw = dofArray(node_left, 3);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(phi0_an - u(rw, itime));

    % reference potential for cathode
    rw = dofArray(node_right, 3);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(phi0_ca - u(rw, itime));

else % for other time steps
    
    F = zeros(ndofs, 1);

    % current for electrode material
    rw = dofArray(node_right, 3);
    F(rw) = -coefs.i_ca * coefs.I_bc;

    Res = Res + F;

    % fixed potential for electrode at anode left
    rw = dofArray(node_left, 3);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(0 - u(rw, itime));

end