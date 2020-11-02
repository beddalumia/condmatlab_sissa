N = 16;
nbit = 64;
numtype = sprintf('int%d',nbit);

E = cell(1,N);
%% PBC Allowed Momentum Values %%
for k = -N/2+1:N/2  % Defined modulo 2pi/N
% (Momentum Loop: each k-sector is decoupled from all the others)

%% List of Momentum States %%
M_k = countKstates(k,N,numtype)
s_a = -ones(1,M_k);   % Basis-set list preallocation
R_a = -ones(1,M_k);   % Periodicity list preallocation
a = 0;                % State counter reset
for s = 0:2^N-1
    R = checkstate(s,k,N,numtype);
    if R >= 0 
        a = a + 1;
        s_a(a) = s;
        R_a(a) = R;
    end
end

%% Hamiltonian %%
H_k = zeros(M_k,M_k);
E{k+N/2} = 0;
for a = 1:M_k
    for i = 0:(N-1)
        j = mod(i+1,N);
        si = bitget(s_a(a),i+1,numtype);
        sj = bitget(s_a(a),j+1,numtype);
        if si == sj
            H_k(a,a) = H_k(a,a) + 1/4;
        else
            H_k(a,a) = H_k(a,a) - 1/4;
            s_flipado = bitxor(s_a(a),2^i+2^j,numtype);     % spin-flip
            [r,L] = representative(s_flipado,N,numtype);    % Repr. of spin-flip
            b = find(s_a==r);                               % Position of spin-flip
            if isempty(b) == 0
                 H_k(a,b) = H_k(a,b) + 1/2*sqrt(R_a(a)./R_a(b))*exp(1i*2*pi/N*k*L);
            end
        end
    end 
end
[E{k+N/2}] = eig(H_k);
% Printing
m = k
E{k+N/2}
end