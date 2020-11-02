Nx = 4;
Ny = 4;
N = Nx*Ny;
nbit = 32;
numtype = sprintf('int%d',nbit);

imax = 2^(N-1);
s = randi(imax);
fprintf('s = %d\n', s);
Lattice2D(s,Nx,Ny,numtype);
Ts = Ty(s,Nx,Ny,numtype);
counter = 1;
fprintf('T(%d)s = %d\n', counter, Ts)
Lattice2D(Ts,Nx,Ny,numtype)

while Ts ~= s
    Ts = Ty(Ts,Nx,Ny,numtype);
    counter = counter + 1;
    fprintf('T(%d)s = %d\n', counter, Ts)
    Lattice2D(Ts,Nx,Ny,numtype)
end
fprintf('Counter: %d\n', counter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s = 38245;
% fprintf('s = %d\n', s)
% Ts = Ty(s,Nx,Ny,numtype);
% counter = 1;
% fprintf('T(%d)s = %d\n', counter, Ts)
% 
% while Ts ~= s
%     Ts = Ty(Ts,Nx,Ny,numtype);
%     counter = counter + 1;
%     fprintf('T(%d)s = %d\n', counter, Ts)
% end
% fprintf('Counter: %d\n', counter)