%================= ELEMENT STIFFNESS MATRIX for CATHODE================================
%
function [kel, rel] = elstif_cathode(ncoord, nelnodes, lmncoord, coefs, elu, eldu)

    % elem solution
    % [ phi_e; c_e; phi_s; j];
    
    % coefficient 2 for 2D, 3 for 3D
    % 2*Pi*R / Pi*R^2 = 2
    % 4*Pi*R^2 / (4/3*Pi*R^3)
    svr = 3 * coefs.eps_inc_ca / coefs.r_ca; % m^-1;

    % Assemble the element stiffness
    % subscripts: p for phi_e; c for phi_c; s for phi_s; j for j
    kpp = zeros(nelnodes, nelnodes); kpc = zeros(nelnodes, nelnodes);
    kps = zeros(nelnodes, nelnodes); kpj = zeros(nelnodes, nelnodes);

    kcp = zeros(nelnodes, nelnodes); kcc = zeros(nelnodes, nelnodes);
    kcs = zeros(nelnodes, nelnodes); kcj = zeros(nelnodes, nelnodes);

    ksp = zeros(nelnodes, nelnodes); ksc = zeros(nelnodes, nelnodes);
    kss = zeros(nelnodes, nelnodes); ksj = zeros(nelnodes, nelnodes);

    kjp = zeros(nelnodes, nelnodes); kjc = zeros(nelnodes, nelnodes);
    kjs = zeros(nelnodes, nelnodes); kjj = zeros(nelnodes, nelnodes);

    rp = zeros(nelnodes, 1); rc = zeros(nelnodes, 1);
    rs = zeros(nelnodes, 1); rj = zeros(nelnodes, 1);

    elp = elu(1:nelnodes);
    elc = elu(nelnodes + 1:2 * nelnodes);
    els = elu(2 * nelnodes + 1:3 * nelnodes);
    elj = elu(3 * nelnodes + 1:end);

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

        % Compute shape functions && derivatives wrt local coords
        %
        xi = xilist(:, intpt);
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
        ratio_mat = (coefs.eps_mat_ca)^(coefs.alpha);
        ratio_inc = (coefs.eps_inc_ca)^(coefs.alpha);

        % -------------------------------------------------------------------------------
        % rp, kpp, kpc
        % -------------------------------------------------------------------------------
        % % coefficients for electrolyte
        % ka = coefs.kappa / coefs.ka_ref;

        % ka_eff = ka * ratio_mat;
        % kad_eff = 2 * ka_eff * coefs.R * coefs.T / coefs.F * (1 - coefs.tp);

        % rp = rp + ka_eff * dNdx * (dNdx' * elp) * w(intpt) * dt ...
        %     -kad_eff / (N' * elc) * dNdx * (dNdx' * elc) * w(intpt) * dt ...
        %     -coefs.sour_p * svr * N * (N' * elj) * w(intpt) * dt;

        % kpp = kpp + ka_eff * (dNdx * dNdx') * w(intpt) * dt;

        % kpc = kpc - (kad_eff / (N' * elc) * (dNdx * dNdx') ...
        %     -kad_eff / (N' * elc)^2 * dNdx * (dNdx' * elc) * N') * w(intpt) * dt;
            
        % kpj = kpj - coefs.sour_p * svr * (N * N') * w(intpt) * dt;

        [kappa, kaDeriv] = calcKappa_Deriv(N' * elc);
        ka_eff = kappa / coefs.ka_ref * ratio_mat;
        kad_eff = 2 * ka_eff * coefs.R * coefs.T / coefs.F * (1 - coefs.tp);

        rp = rp + ka_eff * dNdx * (dNdx' * elp) * w(intpt) * dt ...
            -kad_eff / (N' * elc) * dNdx * (dNdx' * elc) * w(intpt) * dt ...
            -coefs.sour_p * svr * N * (N' * elj) * w(intpt) * dt;

        kap_eff = kaDeriv / coefs.ka_ref * ratio_mat;
        kadp_eff = 2 * kap_eff * coefs.R * coefs.T / coefs.F * (1 - coefs.tp);


        kpp = kpp + ka_eff * (dNdx * dNdx') * w(intpt) * dt;

        kpc = kpc + kap_eff * dNdx * (dNdx' * elp) * N' * w(intpt) * dt ...
            - (kad_eff / (N' * elc) * (dNdx * dNdx') ...
            + (kadp_eff / (N' * elc) - kad_eff / (N' * elc)^2) * dNdx * (dNdx' * elc) * N') * w(intpt) * dt;

        kpj = kpj - coefs.sour_p * svr * (N * N') * w(intpt) * dt;

        % -------------------------------------------------------------------------------
        % rc, kcc
        % -------------------------------------------------------------------------------

        df_eff = coefs.df * ratio_mat / coefs.df_ref;

        rc = rc + coefs.dtco_mat * coefs.eps_mat_ca * N * (N' * eldc) * w(intpt) * dt ...
            +df_eff * dNdx * (dNdx' * elc) * w(intpt) * dt ...
            -coefs.sour_c * svr * N * (N' * elj) * w(intpt) * dt;

        kcc = kcc + coefs.dtco_mat * coefs.eps_mat_ca * (N * N') * w(intpt) * dt ...
            +df_eff * (dNdx * dNdx') * w(intpt) * dt;

        kcj = kcj - coefs.sour_c * svr * (N * N') * w(intpt) * dt;

        % -------------------------------------------------------------------------------
        %  rs, kss
        % -------------------------------------------------------------------------------
        sig = coefs.sigca / coefs.sigca_ref;
        sig_eff = sig * ratio_inc;

        rs = rs + sig_eff * dNdx * (dNdx' * els) * w(intpt) * dt ...
            +coefs.sour_ca * svr * N * (N' * elj) * w(intpt) * dt;

        kss = kss + sig_eff * (dNdx * dNdx') * w(intpt) * dt;

        ksj = ksj + coefs.sour_ca * svr * (N * N') * w(intpt) * dt;

    end

    kel = [kpp, kpc, kps, kpj; kcp, kcc, kcs, kcj; ...
            ksp, ksc, kss, ksj; kjp, kjc, kjs, kjj];
    rel = [rp; rc; rs; rj];

end
