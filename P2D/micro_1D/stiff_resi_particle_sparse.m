%====================== Assemble the global stiffness and residual =================
%
function [K, Res] = stiff_resi_particle_sparse (ncoord, ndof, nnode, coords,...
                     nelnodes, connect, dofArray, coefs, u, du)
    
    nelem = size(connect, 1);

    Ie = zeros ( (ndof * nelnodes)^2, nelem) ;
    Je = zeros ( (ndof * nelnodes)^2, nelem) ;
    Xe = zeros ( (ndof * nelnodes)^2, nelem) ;
    Res = zeros(ndof * nnode, 1);

    % ---------------------------------------------------------------
    % Assemble the stiffness matrix and residual for SPE
    % ---------------------------------------------------------------
    %
    %   Loop over all the elements

    for lmn = 1:nelem
        %
        %   Extract coords of nodes, DOF for the current element
        lmncoord = coords(connect(lmn, :), :);

        % all element nodes dofs 
        elemDofs = dofArray(connect(lmn, :), 1); 

        elu = u(elemDofs);
        eldu = du(elemDofs);
        %
        %   element stiffness and residual

        [kel, rel] = elstif_particle(ncoord, ndof, nelnodes, lmncoord, coefs, elu, eldu);
        %
        %   Add the current element stiffness to the global stiffness

		Ie(:,lmn) = repmat( elemDofs, nelnodes*ndof, 1 );

		Je(:,lmn) = repelem( elemDofs, nelnodes*ndof );

		Xe(:,lmn) = reshape(kel, [], 1);

        Res(elemDofs) = Res(elemDofs) + rel;

    end

    K = sparse(Ie, Je, Xe, ndof * nnode, ndof * nnode);
    
end
