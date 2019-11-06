%====================== Assemble the global stiffness and residual =================
%
function [K, Res, microsolu] = stiff_residual(ncoord, ~, ndofs, coords, ...
        nelnodes, connect, nodes_anode, nodes_cathode, dofArray, coefs, u, du, micromesh, microsolu)

    conn_anode = connect.anode;
    conn_separator = connect.separator;
    conn_cathode = connect.cathode;

    nelem_anode = size(conn_anode, 1);
    nelem_separator = size(conn_separator, 1);
    nelem_cathode = size(conn_cathode, 1);

    Res = zeros(ndofs, 1);

    % ---------------------------------------------------------------
    % Assemble the stiffness matrix and residual for anode
    % ---------------------------------------------------------------
    ndof = 4;
    Ie_an = zeros((ndof * nelnodes)^2, nelem_anode);
    Je_an = zeros((ndof * nelnodes)^2, nelem_anode);
    Xe_an = zeros((ndof * nelnodes)^2, nelem_anode);

    %   Loop over all the elements
    %
    for lmn = 1:nelem_anode
        %
        %   Extract coords of nodes, DOF for the current element
        lmncoord = coords(conn_anode(lmn, :), :);

        % all element nodes dofs
        % dofArray(connect(lmn, :), 1); dofArray(connect(lmn, :), 2); ...
        elemDofs = dofArray(conn_anode(lmn, :), :);
        rws = reshape(elemDofs, [], 1);

        elu = u(rws);
        eldu = du(rws);

        %   element stiffness and residual
        [kel, rel] = elstif_anode(ncoord, nelnodes, lmncoord, coefs, elu, eldu);
        %
        %   Add the current element stiffness to the global stiffness

        Ie_an(:, lmn) = repmat(rws, nelnodes * ndof, 1);

        Je_an(:, lmn) = repelem(rws, nelnodes * ndof);

        Xe_an(:, lmn) = reshape(kel, [], 1);

        Res(rws) = Res(rws) + rel;

    end

    Ie_an = reshape(Ie_an, [], 1);
    Je_an = reshape(Je_an, [], 1);
    Xe_an = reshape(Xe_an, [], 1);

    % remove zero elements in Xe_an
    indx = find(Xe_an);
    Ie_an = Ie_an(indx);
    Je_an = Je_an(indx);
    Xe_an = Xe_an(indx);

    % ---------------------------------------------------------------
    % Assemble the stiffness matrix and residual for cathode
    % ---------------------------------------------------------------
    ndof = 4;
    Ie_ca = zeros((ndof * nelnodes)^2, nelem_cathode);
    Je_ca = zeros((ndof * nelnodes)^2, nelem_cathode);
    Xe_ca = zeros((ndof * nelnodes)^2, nelem_cathode);

    %   Loop over all the elements
    %
    for lmn = 1:nelem_cathode
        %
        %   Extract coords of nodes, DOF for the current element
        lmncoord = coords(conn_cathode(lmn, :), :);

        % all element nodes dofs
        % dofArray(connect(lmn, :), 1); dofArray(connect(lmn, :), 2); ...
        elemDofs = dofArray(conn_cathode(lmn, :), :);
        rws = reshape(elemDofs, [], 1);

        elu = u(rws);
        eldu = du(rws);

        %   element stiffness and residual
        [kel, rel] = elstif_cathode(ncoord, nelnodes, lmncoord, coefs, elu, eldu);
        %
        %   Add the current element stiffness to the global stiffness

        Ie_ca(:, lmn) = repmat(rws, nelnodes * ndof, 1);

        Je_ca(:, lmn) = repelem(rws, nelnodes * ndof);

        Xe_ca(:, lmn) = reshape(kel, [], 1);

        Res(rws) = Res(rws) + rel;

    end

    Ie_ca = reshape(Ie_ca, [], 1);
    Je_ca = reshape(Je_ca, [], 1);
    Xe_ca = reshape(Xe_ca, [], 1);

    % remove zero elements in Xe_ca
    indx = find(Xe_ca);
    Ie_ca = Ie_ca(indx);
    Je_ca = Je_ca(indx);
    Xe_ca = Xe_ca(indx);

    % ---------------------------------------------------------------
    % Assemble the stiffness matrix and residual for separator
    % ---------------------------------------------------------------
    ndof = 2;
    Ie_se = zeros((ndof * nelnodes)^2, nelem_separator);
    Je_se = zeros((ndof * nelnodes)^2, nelem_separator);
    Xe_se = zeros((ndof * nelnodes)^2, nelem_separator);

    %   Loop over all the elements
    %
    for lmn = 1:nelem_separator
        %
        %   Extract coords of nodes, DOF for the current element
        lmncoord = coords(conn_separator(lmn, :), :);

        % all element nodes dofs
        % dofArray(connect(lmn, :), 1); dofArray(connect(lmn, :), 2);
        % only phi_e and c_e
        elemDofs = dofArray(conn_separator(lmn, :), 1:2);
        rws = reshape(elemDofs, [], 1);

        elu = u(rws);
        eldu = du(rws);

        %   element stiffness and residual
        [kel, rel] = elstif_separator(ncoord, nelnodes, lmncoord, coefs, elu, eldu);
        %
        %   Add the current element stiffness to the global stiffness

        Ie_se(:, lmn) = repmat(rws, nelnodes * ndof, 1);

        Je_se(:, lmn) = repelem(rws, nelnodes * ndof);

        Xe_se(:, lmn) = reshape(kel, [], 1);

        Res(rws) = Res(rws) + rel;

    end

    Ie_se = reshape(Ie_se, [], 1);
    Je_se = reshape(Je_se, [], 1);
    Xe_se = reshape(Xe_se, [], 1);

    % remove zero elements in Xe_se
    indx = find(Xe_se);
    Ie_se = Ie_se(indx);
    Je_se = Je_se(indx);
    Xe_se = Xe_se(indx);

    % ---------------------------------------------------------------
    % Assemble the stiffness matrix for pore wall flux
    % ---------------------------------------------------------------

    if coefs.itime == 1

        % all nodal flux dofs
        dofs_j = dofArray([nodes_anode; nodes_cathode], 4);
        Res(dofs_j) = u(dofs_j);
        Ie_f = dofs_j;
        Je_f = dofs_j;
        Xe_f = ones(length(dofs_j), 1);

    else
        stiffness_resid_flux;
    end

    % ---------------------------------------------------------------
    % Assemble global stiffness matrix
    % ---------------------------------------------------------------

    Ie = [Ie_an; Ie_ca; Ie_se; Ie_f];
    Je = [Je_an; Je_ca; Je_se; Je_f];
    Xe = [Xe_an; Xe_ca; Xe_se; Xe_f];

    K = sparse(Ie, Je, Xe, ndofs, ndofs);

end
