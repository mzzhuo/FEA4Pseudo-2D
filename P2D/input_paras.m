%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------
% physical constants
% ---------------------------------------------------------------
coefs.F = 96485.33289; % Faraday's constant (C/mol)
coefs.R = 8.31447; % Gas constant
coefs.T = 298.15; % Absolute temperature 25 celsius

coefs.phi_ref = coefs.R * coefs.T / coefs.F;

% ---------------------------------------------------------------
% interfacial kinetics
% ---------------------------------------------------------------
% nominal Reaction rates (A/m^2)*(mol^3/mol)^(1+alpha)
% coefs.k_s = [1e-05, 3e-07] / coefs.F;
coefs.k_s = [1e-10, 3e-11]; %m^2.5 mol^-0.5 s^-1

% maximum conc in cathode active material (mol/m^3)
coefs.c_max = [2.4983e+04, 5.1218e+04];

% ---------------------------------------------------------------
% initial conditions
% ---------------------------------------------------------------
c0_el = 1.0e3; % Initial concentration in electrolyte (mol/m^3)
phi0_el = 0; % (V)

c0_an = 19624; %0.5164 * coefs.c_max(1);
% so that eta = 0
phi0_an = calcUoc(c0_an, 1, coefs); % (V)

c0_ca = 20046; %0.6740 * coefs.c_max(2);
% so that eta = 0
phi0_ca = calcUoc(c0_ca, 2, coefs); % (V)

% ---------------------------------------------------------------
% initial conditions
% ---------------------------------------------------------------
coefs.l_ref = 1.0e-6; % Characteristic length = length of separator (m)

coefs.eps_inc_an = 0.6; % Volume fraction of active particles in anode
coefs.eps_mat_an = 0.3; % Volume fraction of electrolyte in anode
coefs.eps_fil_an = 0.1; % Volume fraction of filler in anode

coefs.eps_inc_ca = 0.5; % Volume fraction of active particles in cathode
coefs.eps_mat_ca = 0.3; % Volume fraction of electrolyte in cathode
coefs.eps_fil_ca = 0.2; % Volume fraction of filler in cathode

coefs.eps_mat_se = 1.0; % Volume fraction of electrolyte in separator

% particle radius in anode
coefs.r_an = 10e-6; % m

% particle radius in cathode
coefs.r_ca = 10e-6; % m

% % volume fraction of electrolyte
% coefs.eps_mat = 0.5;

% % volume fraction of solid active
% coefs.eps_inc = 0.5;

coefs.alpha = 1.5;

% ---------------------------------------------------------------
% transport properties for electrolyte
% ---------------------------------------------------------------
% ionic conductivity
coefs.kappa = calcKappa(c0_el);
coefs.ka_ref = coefs.kappa;

% diffusivity
coefs.df = 2.7877e-10; %5.34e-10;
coefs.df_ref = coefs.df;

% transference number
coefs.tp = 0.4;

% coefficient for conc time diff term
coefs.dtco_mat = 1 / dt * coefs.l_ref^2 / coefs.df_ref;

% coefficient for current flux term
coefs.i_e = coefs.l_ref / coefs.ka_ref;
% coefficient for species flux term
coefs.q_e = coefs.l_ref / coefs.df_ref / coefs.F;

% coefficient for current source term
coefs.sour_p = (coefs.l_ref)^2 / coefs.ka_ref * coefs.F;
% coefficient for species source term
coefs.sour_c = (coefs.l_ref)^2 / coefs.df_ref * (1 - coefs.tp);

% ---------------------------------------------------------------
% transport properties for active material in anode
% ---------------------------------------------------------------
% active material conductivity
coefs.sigan = 100; % s/m
coefs.sigan_ref = coefs.sigan;

% diffusivity
coefs.dan = 3.9e-14; % m^2/s
coefs.dan_ref = coefs.dan;

coefs.dtco_an = 1 / dt * coefs.l_ref^2 / coefs.dan_ref;

coefs.i_an = coefs.l_ref / coefs.sigan_ref;
% coefs.q_an = coefs.l_ref / coefs.dan_ref / coefs.F;

coefs.sour_an = (coefs.l_ref)^2 / coefs.sigan_ref * coefs.F;

% ---------------------------------------------------------------
% transport properties for active material in cathode
% ---------------------------------------------------------------
% active material conductivity
coefs.sigca = 10; % s/m
coefs.sigca_ref = coefs.sigca;

% diffusivity
coefs.dca = 1e-13; % m^2/s
coefs.dca_ref = coefs.dca;

coefs.dtco_ca = 1 / dt * coefs.l_ref^2 / coefs.dca_ref;

coefs.i_ca = coefs.l_ref / coefs.sigca_ref;
% coefs.q_ca = coefs.l_ref / coefs.dca_ref / coefs.F;

coefs.sour_ca = (coefs.l_ref)^2 / coefs.sigca_ref * coefs.F;

% ---------------------------------------------------------------
% applied external current
% ---------------------------------------------------------------
coefs.I_bc = 30; %A/m^2
