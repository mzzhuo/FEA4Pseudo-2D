%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% physical constants:
% 
coefs.F = 96485.33289;         % Faraday's constant (C/mol)
coefs.R = 8.31447;             % Gas constant
coefs.T = 298.15;              % Absolute temperature 25 celsius

% coefs.phi_ref = coefs.R * coefs.T / coefs.F;

% -------------------------------------------------
% for tese
% -------------------------------------------------

% R_s is passed from macro scale
% if l_ref != R_s, REMEMBER to change the gmsh radius to R_s/l_ref 
coefs.l_ref = R_s; 

% -------------------------------------------------
% anode
% -------------------------------------------------
df_an = 3.9e-14; % m^2/s
df_an_ref = df_an;
dtco_an = 1 / dt * coefs.l_ref^2 / df_an_ref;
q_an = 4 * pi * R_s^2 / df_an_ref / coefs.l_ref;

% -------------------------------------------------
% cathode
% -------------------------------------------------
df_ca = 1e-13; % m^2/s
df_ca_ref = df_ca;
dtco_ca = 1 / dt * coefs.l_ref^2 / df_ca_ref;
q_ca = 4 * pi * R_s^2 / df_ca_ref / coefs.l_ref;

% -------------------------------------------------
% cathode
% -------------------------------------------------
if elemTag == 1
    coefs.ds = df_an;
    coefs.ds_ref = df_an_ref;
    coefs.dtco_inc = dtco_an;
    coefs.q_bc = q_an;
elseif elemTag == 2
    coefs.ds = df_ca;
    coefs.ds_ref = df_ca_ref;
    coefs.dtco_inc = dtco_ca;
    coefs.q_bc = q_ca;
else
    error('wrong elem type!')
end

% % -------------------------------------------------
% % active material diffusivity
% coefs.ds = 1e-13; % m^2/s
% coefs.ds_ref = coefs.ds;

% coefs.dtco_inc = 1 / dt * coefs.l_ref^2 / coefs.ds_ref;

% coefs.q_ca = 4 * pi * R_s^2 / coefs.ds_ref / coefs.l_ref;

