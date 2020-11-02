% LJ Parameters
r0 = 3.758/2^(1/6); %Å
e = 99.55; %cm^{-1}

% Equilibrium Lattice Parameters
a = 3.6502; %Å
c = sqrt(8/3)*a;

% LJ Derivatives [NUMERICAL]
l = a;
r = linspace(l,l+1/10^4,10^3);
ndV = (LJpot(r(2))-LJpot(r(1)))/(r(2)-r(1));
r = [l-(r(2)-r(1)),r];
nddV = (LJpot(r(3))+LJpot(r(1))-2*LJpot(r(2)))/(r(2)-r(1))^2;

% LJ Derivatives [ANALYTICAL]
l = a;
dV = 4*e*(6*r0^6/l^7 - 12*r0^12/l^13);
ddV = 4*e*(156*r0^12/l^14 - 42*r0^6/l^8);

% Nearest Neighbors Lists
RA = {[a, 0, 0], -[a,0,0],[a/2, sqrt(3)*a/2,0], -[a/2, sqrt(3)*a/2,0], [-a/2,sqrt(3)*a/2,0], -[-a/2,sqrt(3)*a/2,0]};
RAB = {[0,0,0], [a,0,0], [a/2,sqrt(3)/2*a,0], [0,0,sqrt(8/3)*a], [a,0,sqrt(8/3)*a], [a/2, sqrt(3)/2*a, sqrt(8/3)*a]};
RBA = {-[0,0,0], -[a,0,0], -[a/2,sqrt(3)/2*a,0], -[0,0,sqrt(8/3)*a], -[a,0,sqrt(8/3)*a], -[a/2, sqrt(3)/2*a, sqrt(8/3)*a]};
SB = [a/2,a/(2*sqrt(3)),sqrt(2/3)*a];

% High Symmetry k-points
    G = [0,0,0];
    A = [0,0,pi/c];
    H = [4/3*pi/a,0,pi/c];
    L = [pi/a,-pi/(sqrt(3)*a),pi/c];
    M = [pi/a,-pi/(sqrt(3)*a),0];
    K = [4/3*pi/a,0,0];

% High Symmetry Paths
Nfactor = 100;
Npoints = cell(1,7);
    %G-M
    Npoints{1} = round(Nfactor*norm(G-M));
    t = linspace(0,1,Npoints{1});
    Kx{1} = G(1) + (M(1)-G(1))*t;
    Ky{1} = G(2) + (M(2)-G(2))*t;
    Kz{1} = G(3) + (M(3)-G(3))*t;
    %M-K
    Npoints{2} = round(Nfactor*norm(M-K));
    t = linspace(0,1,Npoints{2});
    Kx{2} = M(1) + (K(1)-M(1))*t;
    Ky{2} = M(2) + (K(2)-M(2))*t;
    Kz{2} = M(3) + (K(3)-M(3))*t;
    %K-G
    Npoints{3} = round(Nfactor*norm(K-G));
    t = linspace(0,1,Npoints{3});
    Kx{3} = K(1) + (G(1)-K(1))*t;
    Ky{3} = K(2) + (G(2)-K(2))*t;
    Kz{3} = K(3) + (G(3)-K(3))*t;
    %G-A
    Npoints{4} = round(Nfactor*norm(G-A));
    t = linspace(0,1,Npoints{4});
    Kx{4} = G(1) + (A(1)-G(1))*t;
    Ky{4} = G(2) + (A(2)-G(2))*t;
    Kz{4} = G(3) + (A(3)-G(3))*t;
    %A-L
    Npoints{5} = round(Nfactor*norm(A-L));
    t = linspace(0,1,Npoints{5});
    Kx{5} = A(1) + (L(1)-A(1))*t;
    Ky{5} = A(2) + (L(2)-A(2))*t;
    Kz{5} = A(3) + (L(3)-A(3))*t;
    %L-H
    Npoints{6} = round(Nfactor*norm(L-H));
    t = linspace(0,1,Npoints{6});
    Kx{6} = L(1) + (H(1)-L(1))*t;
    Ky{6} = L(2) + (H(2)-L(2))*t;
    Kz{6} = L(3) + (H(3)-L(3))*t;
    %H-A
    Npoints{7} = round(Nfactor*norm(H-A));
    t = linspace(0,1,Npoints{7});
    Kx{7} = H(1) + (A(1)-H(1))*t;
    Ky{7} = H(2) + (A(2)-H(2))*t;
    Kz{7} = H(3) + (A(3)-H(3))*t;


freqs = cell(1,6);
freqs{1} = zeros(1,sum([Npoints{:}]));
freqs{2} = zeros(1,sum([Npoints{:}]));
freqs{3} = zeros(1,sum([Npoints{:}]));
freqs{4} = zeros(1,sum([Npoints{:}]));
freqs{5} = zeros(1,sum([Npoints{:}]));
freqs{6} = zeros(1,sum([Npoints{:}]));
counter = 0;
for i = 1:7 % Number of HighSym Lines
    kx = Kx{i};
    ky = Ky{i};
    kz = Kz{i};
    for j = 1:Npoints{i} 
        k = [kx(j),ky(j),kz(j)];
        % Dynamical Matrix
        A0 = 2*dV/l;
        B0 = 2*ddV-A0;
        D = zeros(6,6);
        % In-plane Blocks
        for r = 1:length(RA)
            D(1:3,1:3) = D(1:3,1:3) + (sin(1/2*k*RA{r}'))^2*(A0*eye(3) + B0*(RA{r}'*RA{r})./(norm(RA{r}))^2);
            D(1:3,1:3) = D(1:3,1:3) + A0/2*eye(3) + B0/2*((RAB{r}-SB)'*(RAB{r}-SB))./(norm(RAB{r}-SB))^2;
            D(4:6,4:6) = D(4:6,4:6) + (sin(1/2*k*RA{r}'))^2*(A0*eye(3) + B0*(RA{r}'*RA{r})./(norm(RA{r}))^2);
            D(4:6,4:6) = D(4:6,4:6) + A0/2*eye(3) + B0/2*((RAB{r}-SB)'*(RAB{r}-SB))./(norm(RAB{r}-SB))^2;
        end
        % Inter-plane Blocks
        for r = 1:length(RAB)
            D(1:3,4:6) = D(1:3,4:6) - exp(1i*k*RAB{r}')*(A0/2*eye(3) + B0/2*((RAB{r}-SB)'*(RAB{r}-SB))./(norm(RAB{r}-SB))^2);
            D(4:6,1:3) = D(4:6,1:3) - exp(1i*k*RBA{r}')*(A0/2*eye(3) + B0/2*((RBA{r}+SB)'*(RBA{r}+SB))./(norm(RBA{r}+SB))^2);
        end
        % EigenThings
        counter = counter + 1;
        eigs = eig(D)*(12*1.6*10^(-4))/(39.948*1.66*10^(-27));
        freqs{1}(counter) = sqrt(eigs(1))/10^(12);
        freqs{2}(counter) = sqrt(eigs(2))/10^(12);
        freqs{3}(counter) = sqrt(eigs(3))/10^(12);
        freqs{4}(counter) = sqrt(eigs(4))/10^(12);
        freqs{5}(counter) = sqrt(eigs(5))/10^(12);
        freqs{6}(counter) = sqrt(eigs(6))/10^(12);
    end
end
k = 1:sum([Npoints{:}]);
figure("Name",'HCP Dispersion')
plot(k,freqs{1},'-','LineWidth',1.5); hold on
plot(k,freqs{2},'-.','LineWidth',1.5);
plot(k,freqs{3},':','LineWidth',1.5);
plot(k,freqs{4},'-','LineWidth',1.5);
plot(k,freqs{5},'-.','LineWidth',1.5);
plot(k,freqs{6},':','LineWidth',1.5);
Value = 0;
Tick = cell(1,6);
Tick{1} = 0;
for i = 1:7
    Value = Value + Npoints{i};
    Tick{i+1} = Value;
    xline(Value,'k');
end
xlim([0,sum([Npoints{:}])]);
xticks([Tick{:}])
xticklabels({'\Gamma','M','K','\Gamma','A','L','H','A'})
title('Phonon Dispersion for HCP Solid-Argon');
ylabel('THz');