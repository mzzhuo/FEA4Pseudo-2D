%================= ELEMENT STIFFNESS MATRIX for SPE================================
%
function [kel, rel] = elstif_particle(ncoord, ndof, nelnodes, lmncoord, coefs, elu, eldu)

    % Assemble the element stiffness
    kel = zeros(ndof * nelnodes, ndof * nelnodes);
    rel = zeros(ndof * nelnodes, 1);

    % Set up integration points && weights
    npoints = numberofintegrationpoints(ncoord, nelnodes);
    xilist = integrationpoints(ncoord, nelnodes, npoints);
    w = integrationweights(ncoord, nelnodes, npoints);
    %
    %  Loop over the integration points
    %%
    for intpt = 1:npoints

        % Compute shape functions && derivatives wrt local coords
        xi = xilist(:, intpt);

        N = shapefunctions(nelnodes, ncoord, xi);
        dNdxi = shapefunctionderivs(nelnodes, ncoord, xi);

        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord' * dNdxi;

        dt = det(dxdxi);

        % Convert shape function derivatives:derivatives wrt global coords
        % dxidx = inv(dxdxi);
        % dNdx = dNdxi * dxidx;
        dNdx = dNdxi / dxdxi;

        % Compute the element stiffness matrix
        df = coefs.ds / coefs.ds_ref;

        rel = rel + coefs.dtco_inc * N * (N' * eldu) * 4 * pi * xi^2 * w(intpt) * dt ...
            +df * dNdx * (dNdx' * elu) * 4 * pi * xi^2 * w(intpt) * dt;

        kel = kel + coefs.dtco_inc * (N * N') * 4 * pi * xi^2 * w(intpt) * dt ...
            +df * (dNdx * dNdx') * 4 * pi * xi^2 * w(intpt) * dt;

    end

end
