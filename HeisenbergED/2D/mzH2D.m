Nx = 2;
Ny = 10;
N = Nx*Ny;
nbit = 32;
numtype = sprintf('int%d',nbit);

%% Loop on Sz values %%
E = cell(1,N+1);
for mz = -N/2:N/2
    %% List of Sz sector states %%
    M0 = factorial(N)/(factorial(N/2))^2; % Apriori Matrix dimension
    s_a = -1*ones(1,M0);     % Basis-set preallocation
    a = 0;                   % State counter reset
    for s = 0:(2^N-1)        % Span of complete Hilbert Space
        m = 0;               % Magnetization reset
        for i = 1:N          % Span of 2D cluster sites (see Lattice2D for convention)
            spin = bitget(s,i,numtype)-1/2; % spin = \pm 1/2
            m = m + spin;    % Magnetization increment
        end
        if m == mz        % Fixed-Magnetization selection rule
            a = a + 1;    % State counter increment
            s_a(a) = s;   % State saved in basis-set
        end
    end
    M = a;                  % Aposteriori Matrix dimension
    if M ~= M0
        %fprintf('THERE IS SOME ERROR')
    end
fprintf('%d-magnetization-sector Hilbert Space has been built.\n',mz);    
    %% Hamiltonian %%
    H = zeros(M,M);
    for a = 1:M
        for x = 0:Nx-1
            for y = 0:Ny-1
                % Particle we are focusing on:
                i = x+y*Nx;
                si = bitget(s_a(a),i+1,numtype);
                % Horizontal N.N. of i-th particle:
                jRh = mod(x+1,Nx) + y*Nx;
                sjRh = bitget(s_a(a),jRh+1,numtype);
                if si == sjRh
                    H(a,a) = H(a,a) + 1/4;
                else
                    H(a,a) = H(a,a) - 1/4;
                    s_flipado = bitxor(s_a(a),2^i+2^jRh,numtype);
                    b = find(s_a==s_flipado);
                    H(a,b) = H(a,b) + 1/2;
                end
                % Vertical N.N. of i-th particle:
                jUp = x + mod(y+1,Ny)*Nx;
                sjUp = bitget(s_a(a),jUp+1,numtype);
                if si == sjUp
                    H(a,a) = H(a,a) + 1/4;
                else
                    H(a,a) = H(a,a) - 1/4;
                    s_flipado = bitxor(s_a(a),2^i+2^jUp,numtype);
                    b = find(s_a==s_flipado);
                    H(a,b) = H(a,b) + 1/2;
                end
            end
        end
    end
    fprintf('%d-magnetization-sector MATRIX has been built.\n',mz);    
    E{mz+N/2+1} = eig(H);
    fprintf('%d-magnetization-sector has been DIAGONALIZED.\n',mz);    
%Closing mz loop
end