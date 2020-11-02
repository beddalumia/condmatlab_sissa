Nx = 2;
Ny = 10;
N = Nx*Ny;
nbit = 32;
numtype = sprintf('int%d',nbit);

%% Zero-Sz Sector Selection %%
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

fprintf('Zero-Magnetization Hilbert Space has been built.\n')
Levels = [];

%% PBC Allowed Momentum Values %%
E = cell(Nx,Ny); basisDIM = 0;
ERROR = [];
for kx = 0:Nx-1  % Defined modulo 2pi/Nx
    for ky = 0:Ny-1  % Defined modulo 2pi/Ny
    % (Momentum Loop: each (kx,ky)-sector is decoupled from all the others)
fprintf('Starting (%d,%d)-momentum BASIS construction.\n',kx,ky);
        %% List of Momentum States %%
        M_k = countSzK2Dstates(basis,kx,ky,Nx,Ny,numtype); basisDIM = basisDIM + M_k;
        s_k = -ones(1,M_k);   % Basis-set list preallocation
        N_k = -ones(1,M_k);   % Normalization list preallocation
        a = 0;                % State counter reset
        for a_Sz = 1:M
            Ns = check2Dstate(basis(a_Sz),kx,ky,Nx,Ny,numtype);
            if Ns > 0 
                a = a + 1;
                s_k(a) = basis(a_Sz);
                N_k(a) = Ns;
            end
        end
        if a ~= M_k
            fprintf('THERE IS SOME ERROR on k=%d selection\n',k)
        end
fprintf('Starting (%d,%d)-momentum MATRIX construction.\n',kx,ky);        
        %% Hamiltonian %%
        H_k = zeros(M_k,M_k);
        for a = 1:M_k
            for x = 0:Nx-1
                for y = 0:Ny-1
                    % Particle we are focusing on:
                    i = x+y*Nx;
                    si = bitget(s_k(a),i+1,numtype);
                    % Horizontal N.N. of i-th particle:
                    jRh = mod(x+1,Nx) + y*Nx;
                    sjRh = bitget(s_k(a),jRh+1,numtype);
                    if si == sjRh
                        H_k(a,a) = H_k(a,a) + 1/4;
                    else
                        H_k(a,a) = H_k(a,a) - 1/4;
                        s_flip = bitxor(s_k(a),2^i+2^jRh,numtype);            % spin-flip
                        [r,Lx,Ly] = representative2D(s_flip,Nx,Ny,numtype);   % Repr. of spin-flip
                        b = find(s_k==r);                                     % Position of spin-flip
                        if isempty(b) == 0
                            fx = 2*pi/Nx*kx*Lx; fy = 2*pi/Ny*ky*Ly;
                            H_k(a,b) = H_k(a,b) + 1/2*sqrt(N_k(b)./N_k(a))*exp(-1i*(fx+fy));
                        end
                    end
                    % Vertical N.N. of i-th particle:
                    jUp = x + mod(y+1,Ny)*Nx;
                    sjUp = bitget(s_k(a),jUp+1,numtype);
                    if si == sjUp
                        H_k(a,a) = H_k(a,a) + 1/4;
                    else
                        H_k(a,a) = H_k(a,a) - 1/4;
                        s_flip = bitxor(s_k(a),2^i+2^jUp,numtype);            % spin-flip
                        [r,Lx,Ly] = representative2D(s_flip,Nx,Ny,numtype);   % Repr. of spin-flip
                        b = find(s_k==r);                                     % Position of spin-flip
                        if isempty(b) == 0
                            fx = 2*pi/Nx*kx*Lx; fy = 2*pi/Ny*ky*Ly;
                            H_k(a,b) = H_k(a,b) + 1/2*sqrt(N_k(b)./N_k(a))*exp(-1i*(fx+fy));
                        end
                    end
                end
            end
        end
fprintf('Starting (%d,%d)-momentum DIAGONALIZATION.\n',kx,ky);  
        E{kx+1,ky+1} = eig(H_k);       
    end
end