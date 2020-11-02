Nx = 3;
Ny = 2;
N = Nx*Ny;
nbit = 32;
numtype = sprintf('int%d',nbit);

% Zero-Sz Sector Selection 
M0 = factorial(N)/(factorial(N/2))^2; % Apriori Matrix dimension
basis = -1*ones(1,M0);  % Basis-set preallocation
a = 0;                  % State counter reset
for s = 0:(2^N-1)       % Span of complete Hilbert Space
   m = 0;               % Magnetization reset
   for i = 1:N          % Span of lattice sites
       spin = bitget(s,i,numtype)-1/2; % spin = \pm 1/2
       m = m + spin;    % Magnetization increment
   end
       if m == 0        % Zero-Magnetization selection rule
          a = a + 1;    % State counter increment
          basis(a) = s; % State saved in basis-set
       end
end
M = a;                  % Aposteriori Matrix dimension
if M ~= M0
    fprintf('THERE IS SOME ERROR on Sz=0 selection\n')
end
repr = [];
for a = 1:M
            s = basis(a);
            [r,Lx,Ly] = representative2D(s,Nx,Ny,numtype);
            repr = [repr;r];
end

repr = unique(repr);
repr = sort(repr);