% van der Waals radius (experimental)
r0 = 3.758/2^(1/6); %Å
% Binding Energy
e = 99.55; % cm^{-1}

% Symmetry Selection
sym = 'hcp'; %['scl','bcc','fcc','hcp']

% Number of Cells
Nx = 100;
Ny = 50;
Nz = 50;

Nsteps = 11;
LatticeParameters = (1.065+ 0.005*(0:Nsteps-1))*r0; % [INIT: 1.04 scl | 1.21 bcc | 1.52 fcc | 1.065 hcp]
E1 = zeros(1,Nsteps);
E2 = zeros(1,Nsteps);

for i = 1:Nsteps
    a = LatticeParameters(i); 

[r, Nb, Lx, Ly, Lz] = generateLattice(a,Nx,Ny,Nz,sym);
x = r(1,:); y = r(2,:); z = r(3,:);

% Compute cohesion energy
bigLx = Nx*Lx/2;
bigLy = Ny*Ly/2;
bigLz = Nz*Lz/2;
if sym == 'hcp'
bigLx = bigLx*10;
bigLy = bigLy*10;
bigLz = bigLz*10;   
end
N = length(x)/Nb;
%E1(i) = TotalcohesionEnergy(x,y,z,bigLx,bigLy,bigLz)/(N*Nb);
E2(i) = cohesionEnergy(x,y,z);
end

%scatter(LatticeParameters,E1,'+'); hold on
plot(LatticeParameters,E2,'.k','MarkerSize',15); hold on

fitted = fit(LatticeParameters',E2','poly3');
h = plot(fitted,'-b');
set(h, 'LineWidth',1.5)
set(h, 'Color',[0.4660, 0.6740, 0.1880])