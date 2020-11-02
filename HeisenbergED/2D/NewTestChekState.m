Nx = 4;
Ny = 4;
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
for kx = 0%-Nx/2+1:Nx/2  % Defined modulo 2pi/Nx
    for ky = 0%-Ny/2+1:Ny/2  % Defined modulo 2pi/Ny
        %fileGB = fopen(sprintf('buenas_noches_kx%dky%d.txt',kx,ky),'w');
        %for a = 1:M
            s = 4350;% basis(a);
            Ns = check2DstateGB(s,kx,ky,Nx,Ny,numtype);
            
            if Ns == -1
                %fprintf('Not a Representative\n');
                %fprintf('mx = %d, my = %d, Ns = %f\n', kx,ky,Ns);
            elseif Ns == 0
                %fprintf('Phase-Match Failed\n');
                %fprintf('mx = %d, my = %d, Ns = %f\n', kx,ky,Ns);
                %repr = [repr;s];
            else
                %repr = [repr;s];
                
            end
        %end
        %fclose(fileGB);
    end
    
end