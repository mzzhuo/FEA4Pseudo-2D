% read gmsh input file for micro particle
infile = './micro_1D/inputfiles/particle1d';

micromesh.elemType = 'Line2';

[micromesh.nnode, micromesh.coords, micromesh.connect, micromesh.bound] = read_particle1d_mesh(infile, micromesh.elemType);

u = zeros(micromesh.nnode, timesteps);
u(:, 1) = c0_an;
% solu_anode(i) stores nodes_anode(i)
microsolu.anode = cell(length(nodes_anode), 1);
microsolu.anode(:, :) = {u};

u = zeros(micromesh.nnode, timesteps);
u(:, 1) = c0_ca;
% solu_cathode(i) stores nodes_cathode(i)
microsolu.cathode = cell(length(nodes_cathode), 1);
microsolu.cathode(:, :) = {u};

% save('./micro_1D/microsolu', 'solu_anode', 'solu_cathode');

clear u;
