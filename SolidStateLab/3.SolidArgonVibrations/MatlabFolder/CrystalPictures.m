% Lattice Parameter
a =  3.6502; %Å [3.5734 SC | 4.1314 BCC | 5.1628 FCC | 3.6502 HCP]
% vdW Radius
r0 =  3.758/2; %Å

% Symmetry Selection
sym = 'hcp'; %['scl','bcc','fcc','hcp']

% Number of Cells
Nx = 16
Ny = 8
Nz = 8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r, Nb, Lx, Ly, Lz] = generateLattice(a,Nx,Ny,Nz,sym);
x = r(1,:); y = r(2,:); z = r(3,:);
rOLD = (r' - [min(x), min(y), min(z)])';
xOLD = rOLD(1,:); yOLD = rOLD(2,:); zOLD = rOLD(3,:); 

drawcrystal(xOLD,yOLD,zOLD,a,r0,Lx,Ly,Lz,Nb);
print(gcf,'SC_SuperCell','-fillpage','-dpdf','-r720')