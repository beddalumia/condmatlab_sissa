function E = cohesionEnergy(x,y,z)
% Total Energy Subroutine                                          
    N = length(x);
    j = find(x==0&y==0&z==0); 
    E = 0;
    for i = 1:N
        dx = x(j) - x(i); dy = y(j) - y(i); dz = z(j) - z(i);
        rij = norm([dx,dy,dz]);
        if rij > 0
            E = E + 1/2*LJpot(rij);
        end
    end
end