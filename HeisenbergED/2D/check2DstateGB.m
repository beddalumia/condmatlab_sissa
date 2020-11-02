function Ns = check2DstateGB(s,mx,my,Nx,Ny,numtype)
    t = s; Fx = 1; Fy = 1; class = zeros(1,Nx*Ny-1);
    xspotted = 0; yspotted = 0;
    for y = 0:Ny-1
        for x = 0:Nx-1
            i = x + y*Nx;
            class(i+1) = t;
            if t < s % It's sure that s is not the representative:
               Ns = -1; return; % -> check has to fail
            elseif t == s
                if isempty(find(xspotted==x,1)) == 1
                    phaseX = x*mx/Nx;
                    xspotted = [xspotted,x];
                    Fx = Fx + exp(1i*2*pi*phaseX)  
                end
                if isempty(find(yspotted==y,1)) == 1
                    phaseY = y*my/Ny;
                    yspotted = [yspotted,y];
                    Fy = Fy + exp(1i*2*pi*phaseY)
                end
            end
            t = Tx(t,Nx,Ny,numtype);
        end
        t = Ty(t,Nx,Ny,numtype);
    end
    F = Fx*Fy; D = length(unique(class));
    class = class'
    xspotted = xspotted';
    yspotted = yspotted';
    if F ~= 0 % Now the check on phases:
        Ns = D*norm(F)^2;   % Compatible with k 
    else
        Ns = 0;             % Not compatible
    end 
    %fprintf(file,'%d %d %f %f\n',s,D,F,Ns);
end