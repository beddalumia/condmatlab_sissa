function [r, Nb, Lx, Ly, Lz]  = generateLattice(a,Nx,Ny,Nz,sym)

N = Nx*Ny*Nz;

% Simple Cubic Lattice
if sym == 'scl'
    Rx = a;
    Ry = a;
    Rz = a;
    L = 2*a;    % "unit" cell length
    Nb = 1;     % Number of basis atoms
    x = -1*ones(1,Nb*N);
    y = -1*ones(1,Nb*N);
    z = -1*ones(1,Nb*N);
    atomcounter = 0;
    for i = 0:Nx-1
        for j = 0:Ny-1
            for k = 0:Nz-1
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx;
                y(atomcounter) = j*Ry;
                z(atomcounter) = k*Rz;
            end
        end
    end
    Lx = L;
    Ly = L;
    Lz = L;

% Body Centered Cubic    
elseif sym == 'bcc'
    Rx = a;
    Ry = a;
    Rz = a;
    L = 2*a;    % "unit" cell length
    Nb = 2;     % Number of basis atoms
    Delta = [a/2, a/2, a/2]; % Sublattice hopper
    x = -1*ones(1,Nb*N);
    y = -1*ones(1,Nb*N);
    z = -1*ones(1,Nb*N);
    atomcounter = 0;
    for i = 0:Nx-1
        for j = 0:Ny-1
            for k = 0:Nz-1
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx;
                y(atomcounter) = j*Ry;
                z(atomcounter) = k*Rz;
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx + Delta(1);
                y(atomcounter) = j*Ry + Delta(2);
                z(atomcounter) = k*Rz + Delta(3);
            end
        end
    end
    Lx = L;
    Ly = L;
    Lz = L;

% Face Centered Cubic
elseif sym == 'fcc'
    Rx = a;
    Ry = a;
    Rz = a;
    L = 2*a;    % "unit" cell length
    Nb = 4;     % Number of basis atoms
    Delta1 = [a/2, a/2, 0]; % Sublattice hopper
    Delta2 = [0, a/2, a/2]; % Sublattice hopper
    Delta3 = [a/2, 0, a/2]; % Sublattice hopper
    x = -1*ones(1,Nb*N);
    y = -1*ones(1,Nb*N);
    z = -1*ones(1,Nb*N);
    atomcounter = 0;
    for i = 0:Nx-1
        for j = 0:Ny-1
            for k = 0:Nz-1
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx;
                y(atomcounter) = j*Ry;
                z(atomcounter) = k*Rz;
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx + Delta1(1);
                y(atomcounter) = j*Ry + Delta1(2);
                z(atomcounter) = k*Rz + Delta1(3);
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx + Delta2(1);
                y(atomcounter) = j*Ry + Delta2(2);
                z(atomcounter) = k*Rz + Delta2(3);
                atomcounter = atomcounter+1;
                x(atomcounter) = i*Rx + Delta3(1);
                y(atomcounter) = j*Ry + Delta3(2);
                z(atomcounter) = k*Rz + Delta3(3);
            end
        end
    end
    Lx = L;
    Ly = L;
    Lz = L;
    
% Hexagonal Close Packed
elseif sym == 'hcp'
    c = sqrt(8/3)*a;
    a1 = [a, 0, 0];
    a2 = [a/2, sqrt(3)/2*a, 0];
    a3 = [0, 0, c]; 
    Lx = 2*sqrt(3)*a;           % "unit" cell lengths
    Ly = 2*2*a;                 % "unit" cell lengths
    Lz = 2*c;                   % "unit" cell lengths
    Nb = 2;     % Number of basis atoms
    Delta = 1/3*(a1+a2)+1/2*(a3); % Sublattice hopper
    Nx = 2*Nx+1;
    Ny = 2*Ny+1;
    N = Nx*Ny*Nz;
    x = -1*ones(1,Nb*N);
    y = -1*ones(1,Nb*N);
    z = -1*ones(1,Nb*N);
    atomcounter = 0;
    for i = 0:Nx-1
        for j = 0:Ny-1
            for k = 0:Nz-1
                atomcounter = atomcounter+1;
                x(atomcounter) = i*a1(1) + j*a2(1) + k*a3(1);
                y(atomcounter) = i*a1(2) + j*a2(2) + k*a3(2);
                z(atomcounter) = i*a1(3) + j*a2(3) + k*a3(3);
                atomcounter = atomcounter+1;
                x(atomcounter) = i*a1(1) + j*a2(1) + k*a3(1) + Delta(1);
                y(atomcounter) = i*a1(2) + j*a2(2) + k*a3(2) + Delta(2);
                z(atomcounter) = i*a1(3) + j*a2(3) + k*a3(3) + Delta(3);
            end
        end
    end
    
% Error Message    
else 
    fprintf('ERROR: Bravais Structure not recognized!\n');
end

% Embeed all coordinates
r = [x; y; z]; % i.e. r(1,:) == x and same for y <-> 2 and z <-> 3

end