function V = LJpot(r)
% Lennard-Jones Potential:

    % van der Waals radius (experimental)
    r0 = 3.758/2^(1/6); %Å
    % Binding Energy
    e = 99.55; %cm^{-1}
    Repulsive = 4*e.*(r0./r).^12;
    Attractive = -4*e.*(r0./r).^6;
    V = Repulsive + Attractive;
    
end