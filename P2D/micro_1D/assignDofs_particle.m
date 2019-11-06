%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% mesh parameters
% 1D problem, in the radian direction of particle
ncoord = 1;

if strcmp(elemType, 'Line2')
    nelnodes = 2;
    nfacenodes = 1;
elseif strcmp(elemType, 'Line3')
    nelnodes = 3;
    nfacenodes = 1;
end

ndof = 1; % each node has 1 dof, concentration

% Assign dofArray: An array of dof indices: c
dofArray = zeros(nnode, ndof);

% dofs for c
dofArray(1:nnode, 1) = 1:nnode;

% total dofs
ndofs = ndof * nnode;
