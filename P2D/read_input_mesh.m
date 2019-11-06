%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2019, Delft), m.zhuo@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==============================Function read_input_file =================
function [nnode, coords, conn, bound] = read_input_mesh_1d(infile, elemType)
    %
    % read gmsh input file
    %
    infile = [infile, '.msh'];
    infile = fopen(infile, 'r');

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

    % ncoord = 1 is 1D problem, 2 for 2D
    coords = coords(:, 1:ncoord);
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

    % get boundary elements
    if ( string(phyname{1,2}) == '"left"' )
        bound.leftElems = connect((connect(:,1) == phyname{1,1}),3:2+nfacenodes);
    else
        fprintf('Left boundary physical number not found!\n');
    end 
    if ( string(phyname{2,2}) == '"right"' )
        bound.rightElems = connect((connect(:,1) == phyname{2,1}),3:2+nfacenodes);
    else
        fprintf('Right boundary physical number not found!\n');
    end 

    % Get nodes of domains:  electrodes and separator
    if (string(phyname{3, 2}) == '"anode"')
        conn.anode = connect((connect(:, 1) == phyname{3, 1}), 3:2 + nelnodes); % SPE
    else
        fprintf('anode domain physical number not found!\n');
    end

    if (string(phyname{4, 2}) == '"separator"')
        conn.separator = connect((connect(:, 1) == phyname{4, 1}), 3:2 + nelnodes); % SPE
    else
        fprintf('separator domain physical number not found!\n');
    end

    if (string(phyname{5, 2}) == '"cathode"')
        conn.cathode = connect((connect(:, 1) == phyname{5, 1}), 3:2 + nelnodes); % SPE
    else
        fprintf('cathode domain physical number not found!\n');
    end

    % get boundary nodes

    bound.left = unique(bound.leftElems);
    bound.right = unique(bound.rightElems);

    % find anode right, cathode left
    nodes_anode = unique(conn.anode);
    nodes_cathode = unique(conn.cathode);

    x_min = min(coords(:, 1));
    x_max = max(coords(:, 1));
    lx_rve = x_max - x_min;
    epsilon1 = lx_rve / 1e8;
    
    maxcoor = max(coords(nodes_anode, 1));
    bound.anright = find(abs(coords(:, 1) - maxcoor) < epsilon1);

    mincoor = min(coords(nodes_cathode, 1));
    bound.caleft = find(abs(coords(:, 1) - mincoor) < epsilon1);

    fclose(infile);
end

