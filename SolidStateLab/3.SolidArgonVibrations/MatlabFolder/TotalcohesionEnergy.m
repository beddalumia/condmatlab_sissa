function E = TotalcohesionEnergy(x,y,z,Lx,Ly,Lz)
% Total Energy Subroutine
    N = length(x);
    E = 0;
    for i = 1:N
        for j = (i+1):N
            dx = x(j) - x(i); dy = y(j) - y(i); dz = z(j) - z(i);
            dx = dx-Lx*round(dx/Lx); dy = dy-Ly*round(dy/Ly); dz = dz-Lz*round(dz/Lz); %(PBCs...)%
            rij = norm([dx,dy,dz]);
            if i~=j
            E = E + LJpot(rij);
            end
        end
    end
end