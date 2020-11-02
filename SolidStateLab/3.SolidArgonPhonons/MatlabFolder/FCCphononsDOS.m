% LJ Parameters
r0 = 3.758/2^(1/6); %Å
e = 99.55; %cm^{-1}

% Equilibrium Lattice Parameter
a = 5.75 %5.1628; %Å

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

% Cubic Unit-Cell in Reciprocal Space
b1 = 2*pi/a*[1, -1, 1]; % 
b2 = 2*pi/a*[1, 1, -1]; % \vec{k} = k1\vec{b1} + k2\vec{b2} + k3\vec{b3}
b3 = 2*pi/a*[-1, 1, 1]; %
Npoints = 250;
k1 = linspace(0,1,Npoints);
k2 = linspace(0,1,Npoints);
k3 = linspace(0,1,Npoints);
% Translating Back to XYZ frame
kpoint = cell(1,Npoints^3);
kCounter = 0;
for i = 1:Npoints
    for j=1:Npoints
        for l=1:Npoints
            kCounter = kCounter + 1;
            kpoint{kCounter} = k1(i).*b1 + k2(j).*b2 + k3(l).*b3;
            %scatter3(kpoint{kCounter}(1),kpoint{kCounter}(2),kpoint{kCounter}(3)); hold on
        end
    end
end

freqs = cell(1,3);
freqs{1} = zeros(1,Npoints^3);
freqs{2} = zeros(1,Npoints^3);
freqs{3} = zeros(1,Npoints^3);
counter = 0;
for j = 1:Npoints^3 
    k = kpoint{j};
    % Dynamical Matrix
    A = 2*dV/a;
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
 
% Density of Modes
figure("Name", 'DOS')
number_of_bins = 50;
bin = linspace(min([freqs{:}]),max([freqs{:}]),number_of_bins+1);
DOS = cell(1,3);
edges = cell(1,3);
for i = 1:3
[DOS{i},edges{i}] = histcounts(freqs{i},bin);
    C = zeros(1,number_of_bins);
    for j = 1:number_of_bins 
        C(j) = (edges{i}(j)+edges{i}(j+1))/2; % Bin centers
    end
    DOS{i} = DOS{i}./(length(freqs{i}))
    %plot(C,DOS{i},'LineWidth',1.5); hold on
end
TOTALDOS = DOS{1}+DOS{2}+DOS{3};
plot(C,TOTALDOS,'LineWidth',2);