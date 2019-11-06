function [c_ss, dcssdJ, u_up] = micro_particle(itime, dt, pwflux, R_s, micromesh, u, elemTag)

    % c_ss = 0.5164 * 2.4984e+04;
    % dcssdJ = 0;
    nnode = micromesh.nnode;
    coords = micromesh.coords;
    connect = micromesh.connect;
    bound = micromesh.bound;
    elemType = micromesh.elemType;

    % u = solu_anode{inode};

    assignDofs_particle;

    input_paras_particle;

    % iteration steps and tolerance
    NR.iter = 10;
    NR.tol = 1e-12;

    % Initialize solution at time step itime with previous step value
    if itime > 1
        u(:, itime) = u(:, itime - 1);
    end

    % solution(concentration) increament of current time step
    % difference bwt time steps: itime - (itime - 1)
    du = zeros(ndofs, 1);

    % Newton-Raphson iterations
    nit = 0;
    relCor = 1.0;

    while ((relCor > NR.tol) && (nit < NR.iter))% Newton Raphson loop
        %
        nit = nit + 1; % iteration number +1

        % Initialize Global Tangent Matrix K and Global RHS
        [K, Res] = stiff_resi_particle_sparse(ncoord, ndof, nnode, coords, nelnodes, ...
            connect, dofArray, coefs, u(:, itime), du);

        applyBCs_particle;

        % K = full(K);

        % solve for the increment
        correc = -K \ Res;

        du = du + correc;
        % update solution at current time step
        u(:, itime) = u(:, itime) + correc;

        % check convergence
        if nit == 1
            correc1 = correc;
        end

        relCor = norm(correc) / norm(correc1);
        % resd = norm(Res) / ndofs;
        % fprintf('ite num %d, rel correc %e, residual %e\n', nit, relCor, resd);

    end

    % solu_anode(inode) = {u};
    u_up = u;

    calcDerivative;

end
