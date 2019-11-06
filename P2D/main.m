%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: fem to solve P2D model
%
% By: Mingzhao ZHUO (Sep. 2019, Delft), mzzhuo@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
tic
clear;
clc;
% close all

% microscale simulation codes are stored there
addpath('./utilities', './micro_1D');

% read gmsh input file
% 1D
% infile = './inputfiles/macro1d';
% elemType = 'Line2';

% 2D triangle
% infile = './inputfiles/macro2d';
% elemType = 'Tri3';

% 2D quad
infile = './inputfiles/macro2d';
elemType = 'Quad4';

[nnode, coords, connect, bound] = read_input_mesh(infile, elemType);

%% Initialization
% initial start time
t_sta = 0.0;

% time step
dt = 1.0;

% store and pass the timestep to micro problem
coefs.dt = dt;

% Number of time steps
% 0, dt, 2*dt, ... n*dt.
n_time = 10; % total time = n_time * dt
timesteps = n_time + 1; % one more to consider the first step

% Assign dofs for macroscale mesh
assignDofs;

% -------------------------------------------------------------------------------
% input material parameters for macro problem
% -------------------------------------------------------------------------------
input_paras;

% -------------------------------------------------------------------------------
% initiate setting of micro particle mesh/solution
% -------------------------------------------------------------------------------
initiateMicro;

% iteration steps and tolerance
NR.iter = 15;
NR.tol = 1e-12;

u = zeros(ndofs, timesteps);
% initial concentration in electrolyte
u(dofsC_e, 1) = c0_el;

% initial guess for iteration, proper guess is important
u(dofsPhi_e, 1) = phi0_el;
u(dofsPhi_an, 1) = phi0_an;
u(dofsPhi_ca, 1) = phi0_ca;

% u_pre = u(:, 1);

% time marching
for itime = 1:timesteps

    % pass the current time step to micro
    coefs.itime = itime;

    % Initialize solution at time step itime with previous step value
    if itime > 1
        u(:, itime) = u(:, itime - 1);
        % u_pre = u(:, itime - 1);
    end

    %
    % solution (concentration) increament of current time step
    % difference bwt time steps: itime - (itime - 1)
    du = zeros(ndofs, 1);

    % Newton-Raphson iterations
    corr = 1.;
    nit = 0;
    resd = 1.0;

    while ((resd > NR.tol) && (nit < NR.iter))% Newton Raphson loop

        nit = nit + 1; % iteration number +1

        % Initialize Global Tangent Matrix K and Global RHS
        [K, Res, microsolu] = stiff_residual(ncoord, nnode, ndofs, coords, nelnodes, ...
            connect, nodes_anode, nodes_cathode, dofArray, coefs, u(:, itime), du, micromesh, microsolu);

        % apply boundary conditions
        applyBCs;

        % solve for the increament
        correc = -K \ Res;

        % update difference between two timesteps
        du = du + correc;
        % update solution at current time step
        u(:, itime) = u(:, itime) + correc;

        % check convergence
        if nit == 1
            correc1 = correc;
        end

        normcorr1 = max(norm(correc1), 1e-20);
        corr = norm(correc) / normcorr1;
        resd = norm(Res) / ndofs;

        fprintf('step %d/%d, ite num %d, corr %e, resd %e\n', ...
            itime, timesteps, nit, corr, resd);
    end

end

% calculate surface concentration c_ss(x, t), overpotential eta(x, t) etc.
calcSecVars;

% load('./c_ss.mat')
% load('./micro_1D/microsolu.mat')

toc

%% output to paraview
% this only for 2D
% Phi = u(dofArray(1:nnode, 1), end);
% C = u(dofArray(1:nnode, 2), end);
% vtfile = ['./results/', 'fieldprofile', '_', 'end'];
% allElems = [connect.anode; connect.separator; connect.cathode];
% toParaview(coords, allElems, elemType, vtfile, Phi, C);
