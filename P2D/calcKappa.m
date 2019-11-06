function kappa = calcKappa(c)

    % coeffieients from Zhang, Du
    a0 = 0.0911;
    a1 = 1.9101e-3;
    a2 = -1.052e-6;
    a3 = 0.1554e-9;

    kappa      = a0 + a1 * c + a2 * c.^2 + a3 * c.^3;

end


