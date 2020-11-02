function [r,Lx,Ly] = representative2D(s,Nx,Ny,numtype)
    r = s;  %
    t = s;  % Inizializations
    Lx = 0; % ---------------
    Ly = 0; % 
    for y = 0:Ny-1
        for x = 0:Nx-1
            if t < r
                r = t;
                Lx = x;
                Ly = y;
            end
            t = Tx(t,Nx,Ny,numtype);
        end
        t = Ty(t,Nx,Ny,numtype);
    end
end