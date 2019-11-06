% -------------------------------------------------------------------------------
% calculate surface concentration c_ss(x, t) and average concentration
% in particle c_av(x, t)
% -------------------------------------------------------------------------------
dofs_bound = micromesh.bound.right;

c_ss_an = zeros(length(nodes_anode), timesteps);
c_av_an = zeros(length(nodes_anode), timesteps);
c_ss_ca = zeros(length(nodes_cathode), timesteps);
c_av_ca = zeros(length(nodes_cathode), timesteps);

for itime = 1:timesteps

    for i = 1:length(nodes_anode)
        % pass macro node, current time step and timestep, pore wall flux
        u_micro = microsolu.anode{i};
        c_ss_an(i, itime) = u_micro(dofs_bound, itime);
        c_av_an(i, itime) = mean(u_micro(:, itime));
    end

    for i = 1:length(nodes_cathode)
        % pass macro node, current time step and timestep, pore wall flux
        u_micro = microsolu.cathode{i};
        c_ss_ca(i, itime) = u_micro(dofs_bound, itime);
        c_av_ca(i, itime) = mean(u_micro(:, itime));
    end

end

% -------------------------------------------------------------------------------
% calculate over potential eta
% -------------------------------------------------------------------------------
% for each node in anode, cathode, assign another dof to eta

% already defined
% dofsPhi_an = dofArray(nodes_anode, 3);
% dofsPhi_ca = dofArray(nodes_cathode, 3);

dofsPhi_e_an = dofArray(nodes_anode, 1);
dofsPhi_e_ca = dofArray(nodes_cathode, 1);

eta_an = u(dofsPhi_an, :) - u(dofsPhi_e_an, :) - calcUoc(c_ss_an, 1, coefs);
eta_ca = u(dofsPhi_ca, :) - u(dofsPhi_e_ca, :) - calcUoc(c_ss_ca, 2, coefs);

% -------------------------------------------------------------------------------
% calculate i0
% -------------------------------------------------------------------------------

dofsC_e_an = dofArray(nodes_anode, 2);
dofsC_e_ca = dofArray(nodes_cathode, 2);
c_e_an = u(dofsC_e_an, :);
c_e_ca = u(dofsC_e_an, :);

i0_an = calcJ0(c_e_an, c_ss_an, 1, coefs) * coefs.F;
i0_ca = calcJ0(c_e_ca, c_ss_ca, 1, coefs) * coefs.F;

% -------------------------------------------------------------------------------
% calculate local current flux
% -------------------------------------------------------------------------------

dofsFlux_an = dofArray(nodes_anode, 4);
dofsFlux_ca = dofArray(nodes_cathode, 4);
curr_an = u(dofsFlux_an, :) * coefs.F;
curr_ca = u(dofsFlux_ca, :) * coefs.F;

% -------------------------------------------------------------------------------
% check
% -------------------------------------------------------------------------------
% curr_an - i0_an .* calcBV(eta_an);