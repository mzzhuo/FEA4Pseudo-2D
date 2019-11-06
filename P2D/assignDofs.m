%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% mesh parameters
if strcmp(elemType, 'Line2')
    nelnodes = 2; 
    ncoord = 1;
    nfacenodes = 1; 
elseif strcmp(elemType, 'Line3')
    nelnodes = 3; 
    ncoord = 1;
    nfacenodes = 1; 
elseif strcmp(elemType, 'Tri3')
    nelnodes = 3; 
    ncoord = 2;
    nfacenodes = 2; 
elseif strcmp(elemType, 'Quad4')
    nelnodes = 4; 
    ncoord = 2;
    nfacenodes = 2; 
end

% each node has 4 dofs:
% potential and concentration of electrolyte
% potential in solid
% pore wall flux mol / (tm^2)
ndof = 4;

% Assign dofArray: An array of dof indices: phi c
% number of lagrange multipliers = no. of nodes
% ndofs = ndof * nnode + nnode;
% nlam = nnode;

% dof_sta = ndof * nnode;

% ndofs = ndof * nnode;

% dofArray =
%  phi_e    c_e      phi_s      j
% [ 1       nnode+1  2*nnode+1  3*nnode+1
%   2       nnode+2  2*nnode+2
%   ...     ...
%   nnode   2*nnode  3*nnode
% ]
%
% dofArray = reshape((1:ndof * nnode), nnode, ndof);

nodes_anode = unique(connect.anode);
nodes_separator = unique(connect.separator);
nodes_cathode = unique(connect.cathode);

dofArray = zeros(nnode, ndof);

% for phi_e at all nodes
dofArray(:, 1) = 1:nnode;

% for phi_c at all nodes
dofArray(:, 2) = nnode + (1:nnode);

% for phi_s at nodes in anode and cathode,
% including interface nodes shared with electrolyte
dofArray(nodes_anode, 3) = 2 * nnode + (1:length(nodes_anode));
ndofs = 2 * nnode + length(nodes_anode);

dofArray(nodes_cathode, 3) = ndofs + (1:length(nodes_cathode));
ndofs = ndofs + length(nodes_cathode);

% for j at nodes in anode and cathode,
% including interface nodes shared with electrolyte
dofArray(nodes_anode, 4) = ndofs + (1:length(nodes_anode));
ndofs = ndofs + length(nodes_anode);

dofArray(nodes_cathode, 4) = ndofs + (1:length(nodes_cathode));
ndofs = ndofs + length(nodes_cathode);

% potential dofs electrolyte
dofsPhi_e = dofArray(1:nnode, 1);

% concentration dofs electrolyte
dofsC_e = dofArray(1:nnode, 2);

% % potential dofs solid particle
dofsPhi_an = dofArray(nodes_anode, 3);
dofsPhi_ca = dofArray(nodes_cathode, 3);

dofsFlux = dofArray([nodes_anode; nodes_cathode], 4);
