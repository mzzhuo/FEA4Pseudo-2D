
% load('./micro_1D/microsolu.mat', 'solu_anode','solu_cathode');
% load('./c_ss', 'c_ss_an', 'c_ss_ca');
% ---------------------------------------------------------------
% for anode
% ---------------------------------------------------------------
infc = 1;

nnode_an = length(nodes_anode);

dof_phie = dofArray(nodes_anode, 1);
dof_ce = dofArray(nodes_anode, 2);
dof_phis = dofArray(nodes_anode, 3);
dof_j = dofArray(nodes_anode, 4);

phi_e = u(dof_phie);
c_e = u(dof_ce);
phi_s = u(dof_phis);
nodeflux = u(dof_j);

c_ss = zeros(nnode_an, 1);

dcssdJ = zeros(nnode_an, 1);

% use anode material properties

elemTag = 1;

for i = 1:nnode_an
    % pass macro node, current time step and timestep, pore wall flux
    u_micro = microsolu.anode{i};
    [a, b, u_up] = micro_particle(coefs.itime, coefs.dt, nodeflux(i), coefs.r_an, micromesh, u_micro, elemTag);
    c_ss(i) = a;
    dcssdJ(i) = b;
    microsolu.anode(i) = {u_up};
    % store c_ss_an
    % c_ss_an(i,coefs.itime) = a;
end

j0 = calcJ0(c_e, c_ss, infc, coefs);
Uoc = calcUoc(c_ss, infc, coefs);
eta = phi_s - phi_e - Uoc;
BV = calcBV(eta);

% Res(dof_j) =  Res(dof_j) + nodeflux - j0 * BV;
Res(dof_j) = nodeflux - j0 .* BV;

[j0wrtc_e, j0wrtc_ss] = calcJ0derivatives(c_e, c_ss, infc, coefs);

bvderiv = calcBVderivatives(eta);

dUocdc_ss = calcUocderivatives(c_ss, infc, coefs);

Ie_fan = repmat(dof_j, 1, 4);
Je_fan = [dof_phie, dof_ce, dof_phis, dof_j];

Xe_f_phie = -j0 .* bvderiv * (-1);
Xe_f_ce = -j0wrtc_e .* BV;
Xe_f_phis = -j0 .* bvderiv;
Xe_f_j = ones(nnode_an, 1) - (j0wrtc_ss .* dcssdJ .* BV + j0 .* bvderiv .* (-dUocdc_ss .* dcssdJ));

Xe_fan = [Xe_f_phie, Xe_f_ce, Xe_f_phis, Xe_f_j];

% ---------------------------------------------------------------
% for cathode
% ---------------------------------------------------------------
infc = 2;

nnode_ca = length(nodes_cathode);

dof_phie = dofArray(nodes_cathode, 1);
dof_ce = dofArray(nodes_cathode, 2);
dof_phis = dofArray(nodes_cathode, 3);
dof_j = dofArray(nodes_cathode, 4);

phi_e = u(dof_phie);
c_e = u(dof_ce);
phi_s = u(dof_phis);
nodeflux = u(dof_j);

c_ss = zeros(nnode_ca, 1);
dcssdJ = zeros(nnode_ca, 1);

% use cathode material properties
elemTag = 2;

for i = 1:nnode_ca
    % pass macro node, current time step and timestep, pore wall flux
    u_micro = microsolu.cathode{i};
    [a, b, u_up] = micro_particle(coefs.itime, coefs.dt, nodeflux(i), coefs.r_ca, micromesh, u_micro, elemTag);
    c_ss(i) = a;
    dcssdJ(i) = b;
    microsolu.cathode(i) = {u_up};
    % store c_ss_ca
    % c_ss_ca(i,coefs.itime) = a;
end

j0 = calcJ0(c_e, c_ss, infc, coefs);
Uoc = calcUoc(c_ss, infc, coefs);
eta = phi_s - phi_e - Uoc;
BV = calcBV(eta);

% Res(dof_j) =  Res(dof_j) + nodeflux - j0 * BV;
Res(dof_j) = nodeflux - j0 .* BV;

[j0wrtc_e, j0wrtc_ss] = calcJ0derivatives(c_e, c_ss, infc, coefs);

bvderiv = calcBVderivatives(eta);

dUocdc_ss = calcUocderivatives(c_ss, infc, coefs);

Ie_fca = repmat(dof_j, 1, 4);
Je_fca = [dof_phie, dof_ce, dof_phis, dof_j];

Xe_f_phie = -j0 .* bvderiv * (-1);
Xe_f_ce = -j0wrtc_e .* BV;
Xe_f_phis = -j0 .* bvderiv;
Xe_f_j = ones(nnode_an, 1) - (j0wrtc_ss .* dcssdJ .* BV + j0 .* bvderiv .* (-dUocdc_ss .* dcssdJ));

Xe_fca = [Xe_f_phie, Xe_f_ce, Xe_f_phis, Xe_f_j];

% ---------------------------------------------------------------
% assemble stiffness
% ---------------------------------------------------------------
Ie_fan = reshape(Ie_fan, [], 1);
Je_fan = reshape(Je_fan, [], 1);
Xe_fan = reshape(Xe_fan, [], 1);

Ie_fca = reshape(Ie_fca, [], 1);
Je_fca = reshape(Je_fca, [], 1);
Xe_fca = reshape(Xe_fca, [], 1);

Ie_f = [Ie_fan; Ie_fca];
Je_f = [Je_fan; Je_fca];
Xe_f = [Xe_fan; Xe_fca];
