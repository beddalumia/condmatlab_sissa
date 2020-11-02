% LJ Parameters
r0 = 3.758/2^(1/6); %Å
e = 99.55; %cm^{-1}

% Equilibrium Lattice Parameter
a = 5.1628; %Å

% LJ Derivatives [NUMERICAL]
l = a/sqrt(2);
r = linspace(l,l+1/10^4,10^3);
ndV = (LJpot(r(2))-LJpot(r(1)))/(r(2)-r(1));
r = [l-(r(2)-r(1)),r];
nddV = (LJpot(r(3))+LJpot(r(1))-2*LJpot(r(2)))/(r(2)-r(1))^2;

% LJ Derivatives [ANALYTICAL]
l = a/sqrt(2);
dV = 4*e*(6*r0^6/l^7 - 12*r0^12/l^13);
ddV = 4*e*(156*r0^12/l^14 - 42*r0^6/l^8);

% Nearest Neighbors List
R = {[a/2, a/2, 0], [0, a/2, a/2], [a/2, 0, a/2], [-a/2, -a/2, 0], [0, -a/2, -a/2], [-a/2, 0, -a/2], [a/2, -a/2, 0], [0, a/2, -a/2], [a/2, 0, -a/2], [-a/2, a/2, 0], [0, -a/2, a/2], [-a/2, 0, a/2]};

% High Symmetry k-points
    G = [0,0,0];
    X = [2*pi/a,0,0];
    L = [pi/a,pi/a,pi/a];
    W = [2*pi/a,pi/a,0];
    K = [3/2*pi/a,3/2*pi/a,0];

% High Symmetry Paths
Nfactor = 100;
Npoints = cell(1,4);
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
    %W-G
    Npoints{3} = round(Nfactor*norm(W-G));
    t = linspace(0,1,Npoints{3});
    Kx{3} = W(1) + (G(1)-W(1))*t;
    Ky{3} = W(2) + (G(2)-W(2))*t;
    Kz{3} = W(3) + (G(3)-W(3))*t;
    %G-L
    Npoints{4} = round(Nfactor*norm(G-L));
    t = linspace(0,1,Npoints{4});
    Kx{4} = G(1) + (L(1)-G(1))*t;
    Ky{4} = G(2) + (L(2)-G(2))*t;
    Kz{4} = G(3) + (L(3)-G(3))*t;


freqs = cell(1,3);
freqs{1} = zeros(1,sum([Npoints{:}]));
freqs{2} = zeros(1,sum([Npoints{:}]));
freqs{3} = zeros(1,sum([Npoints{:}]));
counter = 0;
for i = 1:4 % Number of HighSym Lines
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
        freqs{1}(counter) = sqrt(eigs(1))/10^(12) * 0.4135665538536/2*pi;
        freqs{2}(counter) = sqrt(eigs(2))/10^(12) * 0.4135665538536/2*pi;
        freqs{3}(counter) = sqrt(eigs(3))/10^(12) * 0.4135665538536/2*pi;
    end
end
k = 1:sum([Npoints{:}]);
figure("Name",'FCC Dispersion')
plot(k,freqs{1},'-','LineWidth',1.5); hold on
plot(k,freqs{2},'-.','LineWidth',1.5);
plot(k,freqs{3},':','LineWidth',1.5);
Value = 0;
Tick = cell(1,6);
Tick{1} = 0;
for i = 1:4
    Value = Value + Npoints{i};
    Tick{i+1} = Value;
    xline(Value,'k');
end
xlim([0,sum([Npoints{:}])]);
xticks([Tick{:}])
xticklabels({'\Gamma','X','W','\Gamma','L'})
title('Phonon Dispersion for FCC Solid-Argon');
ylabel('meV');