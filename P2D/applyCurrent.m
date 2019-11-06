function F = applyCurrent(ncoord, coords, ndofs, nfacenodes, dofArray, coefs, bound)
    %
    % ---------------------------------------------------------------
    %  the current flux in the rigth side of electrode
    % ---------------------------------------------------------------

    F = zeros(ndofs, 1);

    npoints = numberofintegrationpoints(ncoord - 1, nfacenodes);

    xi = zeros(ncoord - 1, 1);
    dxdxi = zeros(ncoord, ncoord - 1);

    xilist = integrationpoints(ncoord - 1, nfacenodes, npoints);
    w = integrationweights(ncoord - 1, nfacenodes, npoints);

    bc = bound.rightElems;

    for iel = 1:size(bc, 1)
        %
        rel = zeros(nfacenodes, 1);

        for intpt = 1:npoints

            xi = xilist(:, intpt);

            %
            N = shapefunctions(nfacenodes, ncoord - 1, xi);
            dNdxi = shapefunctionderivs(nfacenodes, ncoord - 1, xi);

            % Compute the jacobian matrix && its determinant
            lmncoord = coords(bc(iel, :), :);
            dxdxi = lmncoord' * dNdxi;

            if (ncoord == 2)
                dt = sqrt(dxdxi(1, 1)^2 + dxdxi(2, 1)^2);
            elseif (ncoord == 3)
                dt = sqrt(((dxdxi(2, 1) * dxdxi(3, 2)) - (dxdxi(2, 2) * dxdxi(3, 1)))^2 ...
                    +((dxdxi(1, 1) * dxdxi(3, 2)) - (dxdxi(1, 2) * dxdxi(3, 1)))^2 ...
                    +((dxdxi(1, 1) * dxdxi(2, 2)) - (dxdxi(1, 2) * dxdxi(2, 1)))^2);
            end

            for a = 1:nfacenodes
                rel(a) = rel(a) + coefs.I_bc * N(a) * w(intpt) * dt;
            end

        end

        % put elem value to global
        % for phi_s
        rw = dofArray(bc(iel, :), 3);
        F(rw) = F(rw) + rel * coefs.i_ca;

    end

end
