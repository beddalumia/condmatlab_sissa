%% List of lattice parameters of interest %%
a = linspace(5,5.6,250);

%% Cohesion Energy at different lattice parameters %%
CohU = zeros(1,length(a));
for l = 1:length(a)
    % van der Waals radius (experimental)
    r0 = 3.758/2^(1/6); %Å
    % Binding Energy
    e = 99.55; % cm^{-1}
    % Symmetry Selection
    sym = 'fcc';
    % Number of Cells
    Nx = 8;
    Ny = 8;
    Nz = 8;
    Nsteps = 11;
    E = zeros(1,Nsteps);
    for i = 1:Nsteps
        [r, Nb, Lx, Ly, Lz] = generateLattice(a(l),Nx,Ny,Nz,sym);
        x = r(1,:); y = r(2,:); z = r(3,:);
        N = length(x)/Nb;
        CohU(l) = cohesionEnergy(x,y,z);
    end
end

%% Density of States, *intensive*, computed at equilibrium %%
DOS = [4.39040000000000e-05,0.000204928000000000,0.000496640000000000,0.000952448000000000,0.00154572800000000,0.00228160000000000,0.00316736000000000,0.00423910400000000,0.00543219200000000,0.00690790400000000,0.00851648000000000,0.0103412480000000,0.0124872960000000,0.0148213760000000,0.0174940160000000,0.0206319360000000,0.0240764160000000,0.0280175360000000,0.0326915840000000,0.0380090880000000,0.0444478720000000,0.0524551680000000,0.0626758400000000,0.0786915840000000,0.102180224000000,0.105160576000000,0.108192640000000,0.110972032000000,0.114103168000000,0.116946304000000,0.119946880000000,0.122776576000000,0.125741824000000,0.128491648000000,0.110858368000000,0.0751828480000000,0.0747947520000000,0.0740446720000000,0.0727651840000000,0.0705964800000000,0.0667864320000000,0.0610794240000000,0.0498443520000000,0.0729265920000000,0.126208512000000,0.152722688000000,0.163704704000000,0.103349376000000,0.0681634560000000,0.0318310400000000];

%% Frequency Range for Phonons at different lattice parameters %%
maxFreq = zeros(1,length(a));
minFreq = zeros(1,length(a));
for p = 1:length(a)
    % LJ Deriva(p)tives [ANALYTICAL]
    l = a(p)/sqrt(2);
    dV = 4*e*(6*r0^6/l^7 - 12*r0^12/l^13);
    ddV = 4*e*(156*r0^12/l^14 - 42*r0^6/l^8);
    % Nearest Neighbors List
    R = {[a(p)/2, a(p)/2, 0], [0, a(p)/2, a(p)/2], [a(p)/2, 0, a(p)/2], [-a(p)/2, -a(p)/2, 0], [0, -a(p)/2, -a(p)/2], [-a(p)/2, 0, -a(p)/2], [a(p)/2, -a(p)/2, 0], [0, a(p)/2, -a(p)/2], [a(p)/2, 0, -a(p)/2], [-a(p)/2, a(p)/2, 0], [0, -a(p)/2, a(p)/2], [-a(p)/2, 0, a(p)/2]};
    % High Symmetry k-points
        G = [0,0,0];
        X = [2*pi/a(p),0,0];
        L = [pi/a(p),pi/a(p),pi/a(p)];
        W = [2*pi/a(p),pi/a(p),0];
        K = [3/2*pi/a(p),3/2*pi/a(p),0];
    % High Symmetry Paths
    Nfactor = 100;
    Npoints = cell(1,5);
        %G-X
        Npoints{1} = round(Nfactor*norm(G-X));
        t = linspace(0,1,Npoints{1});
        Kx{1} = G(1) + (X(1)-G(1))*t;
        Ky{1} = G(2) + (X(2)-G(2))*t;
        Kz{1} = G(3) + (X(3)-G(3))*t;
        %X-W
        Npoints{2} = round(Nfactor*norm(X-W));
        t = linspace(0,1,Npoints{2});
        Kx{2} = X(1) + (W(1)-X(1))*t;
        Ky{2} = X(2) + (W(2)-X(2))*t;
        Kz{2} = X(3) + (W(3)-X(3))*t;
        %W-K
        Npoints{3} = round(Nfactor*norm(W-K));
        t = linspace(0,1,Npoints{3});
        Kx{3} = W(1) + (K(1)-W(1))*t;
        Ky{3} = W(2) + (K(2)-W(2))*t;
        Kz{3} = W(3) + (K(3)-W(3))*t;
        %K-G
        Npoints{4} = round(Nfactor*norm(K-G));
        t = linspace(0,1,Npoints{4});
        Kx{4} = K(1) + (G(1)-K(1))*t;
        Ky{4} = K(2) + (G(2)-K(2))*t;
        Kz{4} = K(3) + (G(3)-K(3))*t;
        %G-L
        Npoints{5} = round(Nfactor*norm(G-L));
        t = linspace(0,1,Npoints{5});
        Kx{5} = G(1) + (L(1)-G(1))*t;
        Ky{5} = G(2) + (L(2)-G(2))*t;
        Kz{5} = G(3) + (L(3)-G(3))*t;
    freqs = cell(1,3);
    freqs{1} = zeros(1,sum([Npoints{:}]));
    freqs{2} = zeros(1,sum([Npoints{:}]));
    freqs{3} = zeros(1,sum([Npoints{:}]));
    counter = 0;
    for i = 1:5 % Number of HighSym Lines
        kx = Kx{i};
        ky = Ky{i};
        kz = Kz{i};
        for j = 1:Npoints{i} 
            k = [kx(j),ky(j),kz(j)];
            % Dynamical Matrix
            A = 2*dV/l;
            B = 2*ddV-A;
            D = zeros(3,3);
            for r = 1:length(R)
                D = D + ((sin(1/2*k*R{r}'))^2*(A*eye(3) + B*(R{r}'*R{r})./(norm(R{r}))^2));
            end
            % EigenThings
            counter = counter + 1;
            eigs = eig(D)*(12*1.6*10^(-4))/(39.948*1.66*10^(-27));
            freqs{1}(counter) = sqrt(eigs(1))/10^(12);
            freqs{2}(counter) = sqrt(eigs(2))/10^(12);
            freqs{3}(counter) = sqrt(eigs(3))/10^(12);
        end
    end
    maxFreq(p) = max(freqs{3});
    minFreq(p) = min(freqs{1});
end

%% Helmholtz Free Energy Computation and Minimization for several a and T values
F = zeros(30,11);
aEQ = zeros(1,30);
for j = 1:53
    T = (j-1)+0.0001; %Kelvin
    for i = 1:length(a)
        w = linspace(minFreq(i)+0.0001,maxFreq(i),length(DOS));
        % hbar set to one
        h = 1;
        % boltzmann set to one
        k = 1;
        % Helmholtz Free Energy
        F(j,i) = CohU(i) + sum(k*T*log(1-exp(-h*w/(k*T))).*DOS);
    end
    plot(a,F(j,:),'-k'); hold on
    legend off
    aEQ(j) = a(F(j,:) == min(F(j,:)));
    scatter(aEQ(j),min(F(j,:)),5,'x','r','LineWidth',1);
end
xlabel('Lattice Parameter [Å]', 'Interpreter','tex');
ylabel('Helmholtz Free Energy (Classical) [cm^{-1}]', 'Interpreter','tex');
title('Free Energy Minimization for T in range [0, 52] K', 'Interpreter','tex');
figure("Name",'Thermal Expansion')
T = 0:52+0.0001;
fitted = fit(T',aEQ','poly2');
plot(T,aEQ,'+k'); hold on
T = 0:100;
fitY = fitted.p1*T.^2 + fitted.p2*T + fitted.p3;
plot(T,fitY,'--r','LineWidth',1.2);
ylim([5.1, 5.7]);
yline(5.1628,':g','LineWidth', 1.2);
xlabel('Temperature [K]');
ylabel('Equilibrium Lattice Parameter [Å]')
title('Thermal Expansion of FCC Solid-Argon within QHLD');
legend('Computed Data', 'Quadratic Fit', 'Zero-T Value');