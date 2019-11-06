%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==============================Function read_input_file =================
function [nnode, coords, connect, bound] = read_RVE_mesh(infile, elemType)
    %
    % read gmsh input file
    %
    infile = [infile, '.msh'];
    infile = fopen(infile, 'r');

    if strcmp(elemType, 'Line2')
        nelnodes = 2;
    elseif strcmp(elemType, 'Line3')
        nelnodes = 3; % triangle element
    end

    %
    %%
    line = fgetl(infile);
    %
    % get physical line number for interface elems
    while (~strcmp(line, '$PhysicalNames'))
        line = fgetl(infile);
    end

    line = fgetl(infile);
    nphy = sscanf(line, '%d'); % number of physical names
    %
    phyname = cell(nphy, 2);
    %
    for i = 1:nphy
        line = fgetl(infile);
        buf = textscan(line, '%d %d %s');
        phyname{i, 1} = buf{2}; phyname{i, 2} = buf{3};
    end

    %%
    % number of nodes
    while (~strcmp(line, '$Nodes'))
        line = fgetl(infile);
    end

    %
    line = fgetl(infile);
    nnode = sscanf(line, '%d');
    %
    % node coordinates stored in coords
    coords = zeros(nnode, 3); % gmsh 2D/3D mesh both have 3 components

    for i = 1:nnode
        line = fgetl(infile);
        buf = sscanf(line, '%f %f %f %f');
        coords(i, :) = buf(2:4)'; % not store node ID! make sure is 1,2,3....
    end

    coords = coords(:, 1:1); % this is 1D problem
    %%
    %  No. elements && connectivity
    while (~strcmp(line, '$Elements'))
        line = fgetl(infile);
    end

    line = fgetl(infile);
    nelem = sscanf(line, '%d');
    %%
    % gmsh file format:
    % elem ID, elem type, no of tags, tags...,
    % 1 2-node line.
    % 2 3-node triangle.
    % 3 4-node quadrangle.
    %
    connect = zeros(nelem, 5 + nelnodes);

    for i = 1:nelem
        line = fgetl(infile);
        buf = sscanf(line, '%f %f %f %f %f %f %f %f');
        connect(i, 1:length(buf)) = buf; % elem ID not stored!
    end

    %
    connect(:, 1:3) = [];
    %
    %% correct connectivity of elements with negative jacobian
    %%
    % Get nodes of domains: SPEC and electrodes
    %

    if (string(phyname{1, 2}) == '"particle1d"')
        connect = connect((connect(:, 1) == phyname{1, 1}), 3:2 + nelnodes); % SPE
    else
        fprintf('Matrix physical number not found!');
    end

    x_min = min(coords(:, 1));
    x_max = max(coords(:, 1));
    lx_rve = x_max - x_min;

    % get boundary nodes
    epsilon1 = lx_rve / 1e8;

    bound.left = find(abs(coords(:, 1) - x_min) < epsilon1);
    bound.right = find(abs(coords(:, 1) - x_max) < epsilon1);

    fclose(infile);
end
