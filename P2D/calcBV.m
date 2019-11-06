function BV = calcBV( eta )

    F = 96485.33289; % Faraday's constant (C/mol)
    R = 8.31447; % Gas constant
    T = 298.15; % Absolute temperature

    BV = exp(0.5 * F / R / T * eta ) - exp(-0.5 * F / R / T * eta );

end
