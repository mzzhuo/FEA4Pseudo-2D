if nfacenodes == 1
    applyBCs_1d;
elseif nfacenodes == 2
    applyBCs_2d;
else
    error('Wrong BCs!')
end
