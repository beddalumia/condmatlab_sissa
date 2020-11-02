N = 8;
nbit = 16;
numtype = sprintf('int%d',nbit);

%% Zero-Sz Sector Selection %%
M0 = factorial(N)/(factorial(N/2))^2; % Apriori Matrix dimension
basis = -1*ones(1,M0);  % Basis-set preallocation
a = 0;                  % State counter reset
for s = 0:(2^N-1)       % Span of complete Hilbert Space
   m = 0;               % Magnetization reset
   for i = 1:N          % Span of chain sites
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
fprintf('Zero-Magnetization Hilbert Space has been constructed.\n')

%% PBC Allowed Momentum Values %%
E = cell(1,N); H = cell(1,N); Levels = [];
for k = 0:N-1 % Defined modulo 2pi/N
% (Momentum Loop: each k-sector is decoupled from all the others)
fprintf('Starting %d-momentum sector BASIS building.\n',k)
    %% List of Momentum States %%
    M_k = countSzKstates(basis,k,N,numtype);
    s_k = -ones(1,M_k);   % Basis-set list preallocation
    R_k = -ones(1,M_k);   % Periodicity list preallocation
    a = 0;                % State counter reset
    for a_Sz = 1:M
        R = checkstate(basis(a_Sz),k,N,numtype);
        if R >= 0 
            a = a + 1;
            s_k(a) = basis(a_Sz);
            R_k(a) = R;
        end
    end
    if a ~= M_k
        fprintf('THERE IS SOME ERROR on k=%d selection\n',k)
    end
fprintf('Starting %d-momentum sector MATRIX building.\n',k)    
    %% Hamiltonian %%
    H_k = zeros(M_k,M_k);
    for a = 1:M_k
        for i = 0:(N-1)
            j = mod(i+1,N);
            si = bitget(s_k(a),i+1,numtype);
            sj = bitget(s_k(a),j+1,numtype);
            if si == sj
                H_k(a,a) = H_k(a,a) + 1/4;
            else
                H_k(a,a) = H_k(a,a) - 1/4;
                s_fliped = bitxor(s_k(a),2^i+2^j,numtype);     % spin-flip
                [r,L] = representative(s_fliped,N,numtype);    % Repr. of spin-flip
                b = find(s_k==r);                              % Position of spin-flip
                if isempty(b) == 0
                     H_k(a,b) = H_k(a,b) + 1/2*sqrt(R_k(a)./R_k(b))*exp(1i*2*pi/N*k*L);
                end
            end
        end 
    end
fprintf('Starting %d-momentum sector DIAGONALIZATION.\n',k)
    H{k+1} = H_k;
    E{k+1} = eig(H_k);
    Levels = [Levels; sort(real(E{k+1}))];
end