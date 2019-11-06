%%
close all
%% plot the potential in electrolyte
figure
hold on;
plot(coords(:, 1), u(dofsPhi_e, end), 'o', 'MarkerSize', 4);

%% plot the concentration in electrolyte
% figure
hold on;
% plot(coords(:,1),u(dofsC_e,2),'s','MarkerSize',4)
% plot(coords(:,1),u(dofsC_e,9),'o','MarkerSize',4)
plot(coords(:, 1), u(dofsC_e, end), 'o', 'MarkerSize', 4)
axis([0 225 600 1400]);
%%
data = [coords(:, 1), u(dofsC_e, 16)];
data = sortrows(data, 1);
fileID = fopen('fem_30.txt', 'w');
fprintf(fileID, '%12.6e, %12.6e \n', data');
fclose(fileID);

%%
% plot(comsol30(:,1),comsol30(:,2),'-r','Linewidth',2)
plot(comsol1000(:, 1), comsol1000(:, 2), '-r', 'Linewidth', 2)
%%
hold on;
plot(celeft(:, 1), celeft(:, 2), '-r', 'Linewidth', 2)
plot(ceright(:, 1), ceright(:, 2), '-b', 'Linewidth', 2)
%% plot the potential in cathode
figure
hold on;
plot(coords([nodes_anode; nodes_cathode], end), u(dofsFlux, 290) * coefs.F, 's', 'MarkerSize', 4)

% plot(coords(nodes_anode,1),curr_an(:,end),'o','MarkerSize',8)

%% plot the potential in anode
figure
hold on;
plot(coords(nodes_anode, 1), u(dofsPhi_an, end), 's', 'MarkerSize', 4)

%% plot the potential in cathode
figure
hold on;
plot(coords(nodes_cathode, 1), u(dofsPhi_ca, end), 's', 'MarkerSize', 4)

%% plot the potential in cathode
figure
hold on;
plot(coords(nodes_anode, 1), c_ss_an(:, end), 'bo', 'MarkerSize', 4)
plot(coords(nodes_anode, 1), c_av_an(:, end), 'ro', 'MarkerSize', 4)
plot(coords(nodes_anode, 1), c_ss_ca(:, end), 'bs', 'MarkerSize', 4)
plot(coords(nodes_anode, 1), c_av_ca(:, end), 'rs', 'MarkerSize', 4)

%% plot two end conc evolve with time
figure
hold on;
plot(0:dt:n_time * dt, u(dofArray(node_left, 2), :), 'o', 'MarkerSize', 4)
plot(0:dt:n_time * dt, u(dofArray(node_right, 2), :), 'o', 'MarkerSize', 4)

%%
data = [0:dt:n_time * dt, u(dofArray(node_left, 2), :), u(dofArray(node_right, 2), :)];

fileID = fopen('fem_conevol.txt', 'w');
fprintf(fileID, '%12.6e, %12.6e, %12.6e \n', data');
fclose(fileID);

%%
% eta = -0.2:0.01:0.2;
% bv = calcBV(eta);
% figure;
% plot(eta,bv,'o');

%%
% R_s = 1e-6;
% itime = 2;
% u_micro = solu_anode{1};
% pwflux = 0;
% dt = 1;
% [c_ss, dcssdJ, u_up] = micro_particle(itime, dt, pwflux, R_s,micromesh,u_micro, 1);
