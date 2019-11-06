%================= ELEMENT STIFFNESS MATRIX for SEPARATOR================================
%
function [kel, rel] = elstif_separator(ncoord, nelnodes, lmncoord, coefs, elu, eldu)

    % elem solution
    % [ phi_e; c_e; phi_s; j];

    % Assemble the element stiffness
    % subscripts: p for phi_e; c for phi_c; s for phi_s; j for j
    kpp = zeros(nelnodes, nelnodes);
    kpc = zeros(nelnodes, nelnodes);
    kcp = zeros(nelnodes, nelnodes);
    kcc = zeros(nelnodes, nelnodes);

    rp = zeros(nelnodes, 1);
    rc = zeros(nelnodes, 1);

    elp = elu(1:nelnodes);
    elc = elu(nelnodes + 1:2 * nelnodes);

    eldc = eldu(nelnodes + 1:2 * nelnodes);
    %
    %  Set up integration points && weights
    %
    npoints = numberofintegrationpoints(ncoord, nelnodes);
    xilist = integrationpoints(ncoord, nelnodes, npoints);
    w = integrationweights(ncoord, nelnodes, npoints);
    %
    %  Loop over the integration points
    %%
    for intpt = 1:npoints

        %     Compute shape functions && derivatives wrt local coords
        %
        xi = xilist(:, intpt);
        %   for i = 1:ncoord
        %     xi(i) = xilist(i,intpt);
        %   end
        N = shapefunctions(nelnodes, ncoord, xi);
        dNdxi = shapefunctionderivs(nelnodes, ncoord, xi);

        %   Compute the jacobian matrix && its determinant
        %
        dxdxi = lmncoord' * dNdxi;

        dt = det(dxdxi);

        %   Convert shape function derivatives:derivatives wrt global coords
        % dxidx = inv(dxdxi);
        % dNdx = dNdxi * dxidx;
        dNdx = dNdxi / dxdxi;

        % calculate the residual
        ratio_mat = (coefs.eps_mat_se)^(coefs.alpha);

        % -------------------------------------------------------------------------------
        % rp, kpp, kpc
        % -------------------------------------------------------------------------------
        % % coefficients for electrolyte
        % ka = coefs.kappa / coefs.ka_ref;

        % ka_eff = ka * ratio_mat;
        % kad_eff = 2 * ka_eff * coefs.R * coefs.T / coefs.F * (1 - coefs.tp);

        % rp = rp + ka_eff * dNdx * (dNdx' * elp) * w(intpt) * dt ...
        %     -kad_eff / (N' * elc) * dNdx * (dNdx' * elc) * w(intpt) * dt;

        % kpp = kpp + ka_eff * (dNdx * dNdx') * w(intpt) * dt;

        % kpc = kpc - (kad_eff / (N' * elc) * (dNdx * dNdx') ...
        %     -kad_eff / (N' * elc)^2 * dNdx * (dNdx' * elc) * N') * w(intpt) * dt;

        [kappa, kaDeriv] = calcKappa_Deriv(N' * elc);
        ka_eff = kappa / coefs.ka_ref * ratio_mat;
        kad_eff = 2 * ka_eff * coefs.R * coefs.T / coefs.F * (1 - coefs.tp);

        rp = rp + ka_eff * dNdx * (dNdx' * elp) * w(intpt) * dt ...
            -kad_eff / (N' * elc) * dNdx * (dNdx' * elc) * w(intpt) * dt;

        kap_eff = kaDeriv / coefs.ka_ref * ratio_mat;
        kadp_eff = 2 * kap_eff * coefs.R * coefs.T / coefs.F * (1 - coefs.tp);

        kpp = kpp + ka_eff * (dNdx * dNdx') * w(intpt) * dt;

        kpc = kpc + kap_eff * dNdx * (dNdx' * elp) * N' * w(intpt) * dt ...
            -(kad_eff / (N' * elc) * (dNdx * dNdx') ...
            +(kadp_eff / (N' * elc) - kad_eff / (N' * elc)^2) * dNdx * (dNdx' * elc) * N') * w(intpt) * dt;

        % -------------------------------------------------------------------------------
        % rc, kcc
        % -------------------------------------------------------------------------------

        df_eff = coefs.df * ratio_mat / coefs.df_ref;

        rc = rc + coefs.dtco_mat * coefs.eps_mat_se * N * (N' * eldc) * w(intpt) * dt ...
            +df_eff * dNdx * (dNdx' * elc) * w(intpt) * dt;

        kcc = kcc + coefs.dtco_mat * coefs.eps_mat_se * (N * N') * w(intpt) * dt ...
            +df_eff * (dNdx * dNdx') * w(intpt) * dt;

    end

    kel = [kpp, kpc; kcp, kcc];
    rel = [rp; rc];

end
