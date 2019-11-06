function [j0wrtc_e, j0wrtc_ss] = calcJ0derivatives(c_e, c_ss, i, coefs)
    % input c_s, c_e should have unit of A/m^2
    % i = 1 or 2, 
    % 1-anode_spe interface,
    % 2-cathode_spe interface

    k_s = coefs.k_s;

    % maximum conc in cathode active material (mol/m^3)
    c_max = coefs.c_max;

    alpha_a = 0.5;
    alpha_c = 0.5;

    if (i == 1)

        % with respect to c_e, concentration in electrolyte
        j0wrtc_e = k_s(1) * (c_max(1) - c_ss).^alpha_a .* (c_ss).^alpha_c .* alpha_a .* (c_e).^(alpha_a - 1);

        % with respect to c_ss, surface concentration in particle
        j0wrtc_ss = k_s(1) * (c_e).^alpha_a .* (c_ss).^alpha_c .* alpha_a .* (c_max(1) - c_ss).^(alpha_a - 1) * (-1) ...
            +k_s(1) * (c_max(1) - c_ss).^alpha_a .* (c_e).^alpha_a .* alpha_c .* (c_ss).^(alpha_c - 1); % A/m^2


    elseif (i == 2)
        
        % with respect to c_e, concentration in electrolyte
        j0wrtc_e = k_s(2) * (c_max(2) - c_ss).^alpha_a .* (c_ss).^alpha_c .* alpha_a .* (c_e).^(alpha_a - 1);

        % with respect to c_ss, surface concentration in particle
        j0wrtc_ss = k_s(2) * (c_e).^alpha_a .* (c_ss).^alpha_c .* alpha_a .* (c_max(2) - c_ss).^(alpha_a - 1) * (-1) ...
            +k_s(2) * (c_max(2) - c_ss).^alpha_a .* (c_e).^alpha_a .* alpha_c .* (c_ss).^(alpha_c - 1); % A/m^2

    else
        error('wrong input of interface number!');
    end

    %
end
