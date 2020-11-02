% van der Waals radius (experimental)
r0 = 3.758; %Å
% Binding Energy
e = 99.55; % cm^{-1}

% Symmetry Selection
sym = {'scl','bcc','fcc','hcp'};
for i = 1:4

if sym{i} == 'scl'
    L6 = 8.40192397482537;
    L12 = 6.20214904504752;
    gFactor = 1;
elseif sym{i} == 'bcc'
    L6 = 12.2536678672899;
    L12 = 9.11418326807536;
    gFactor = 3^(-1/2)*2;
elseif sym{i} == 'fcc'
    L6 = 14.4539210437416;
    L12 = 12.1318801965446;
    gFactor = 2^(1/2);
elseif sym{i} == 'hcp'
    L6 = 14.4548972778391;
    L12 = 12.1322937690989;
    gFactor = 1;
end

req = (L12/L6)^(1/6)*r0;
Ucoh = 1/2*e*L6^2/L12;
aeq = req * gFactor;

fprintf('\n***************************************\n')
fprintf('For %s lattice:\n',sym{i})
fprintf('Cohesion Energy is: %f\n',Ucoh)
fprintf('Equilibrium Lattice Parameter is: %f\n',aeq)
fprintf('***************************************\n')

end