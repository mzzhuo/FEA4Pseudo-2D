function j0 = calculateJ0(c_e, c_ss, i, coefs)
    % input c_s, c_e should have unit of A/m^2
    % i = 1 or 2,
    % 1-anode_spe interface,
    % 2-cathode_spe interface

    % F = 96485.33289; % Faraday's constant (C/mol)

    % nominal Reaction rates (A/m^2)*(mol^3/mol)^(1+alpha)
    k_s = coefs.k_s;

    % maximum conc in cathode active material (mol/m^3)
    c_max = coefs.c_max;

    alpha_a = 0.5;
    alpha_c = 0.5;

    if (i == 1)

        % error('anode not considered');

        j0 = k_s(1) * (c_max(1) - c_ss).^alpha_a .* (c_e).^alpha_a .* (c_ss).^alpha_c; % A/m^2

        % if ((c_max(1) - c_e) <= 0 || c_e <= 1)
        %     i0 = 0.0001;
        % end

    elseif (i == 2)
        %
        j0 = k_s(2) * (c_max(2) - c_ss).^alpha_a .* (c_e).^alpha_a .* (c_ss).^alpha_c; % A/m^2

    else
        error('wrong input of interface number!');
    end

    %
end
