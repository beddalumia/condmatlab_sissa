% van der Waals radius (experimental)
r0 = 3.758/2^(1/6); %Å
% Binding Energy
e = 99.55; %cm^{-1}
% Lattice Parameter
a =  3.6502; %Å [3.5734 SC | 4.1314 BCC | 5.1628 FCC | 3.6502 HCP]

% Symmetry Selection
sym = 'hcp'; %['scl','bcc','fcc','hcp']

% Number of Cells
E2 = zeros(1,4)
N = E2;
for Nx = 2:2:50
Ny = Nx;
Nz = Nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[rOLD, Nb, Lx, Ly, Lz] = OLDgenerateLattice(a,Nx,Ny,Nz,sym);
[r, Nb, Lx, Ly, Lz] = generateLattice(a,Nx,Ny,Nz,sym);
x = r(1,:); y = r(2,:); z = r(3,:);
%xOLD = rOLD(1,:); yOLD = rOLD(2,:); zOLD = rOLD(3,:); 
% Compute cohesion energy
bigLx = Nx*Lx/2;
bigLy = Ny*Ly/2;
bigLz = Nz*Lz/2;
N(Nx/2) = length(x);
%E1 = TotalcohesionEnergy(xOLD,yOLD,zOLD,bigLx,bigLy,bigLz)/(N*Nb);
E2(Nx/2) = cohesionEnergy(x,y,z);
%scatter(N*Nb,E1,'+','b'); hold on
%scatter(N*Nb,E2,'s','r'); hold on
end
plot(N,E2,'.-'); hold on