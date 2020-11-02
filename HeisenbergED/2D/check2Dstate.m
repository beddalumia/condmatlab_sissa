function Ns = check2Dstate(s,mx,my,Nx,Ny,numtype)
    t = s; class = zeros(1,Nx*Ny-1); F = 0;
    for y = 0:Ny-1
        for x = 0:Nx-1
            i = x + y*Nx;
            class(i+1) = t;
            if t < s % It's sure that s is not the representative:
                Ns = -1; return; % -> check has to fail
            elseif t == s
                phaseX = x*mx/Nx;
                phaseY = y*my/Ny;
                F = F + exp(-1i*2*pi*(phaseX+phaseY));
            end
            t = Tx(t,Nx,Ny,numtype);
        end
        t = Ty(t,Nx,Ny,numtype);
    end
    D = length(unique(class));
    if norm(F) >= 1 % Now the check on phases:
        Ns = D*norm(F)^2;   % Compatible with k
    else
        Ns = 0;             % Not compatible
    end
end