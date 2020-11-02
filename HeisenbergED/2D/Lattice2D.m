function Lattice2D(s,Nx,Ny,numtype)
    fprintf('\n')
    for y = (Ny-1):-1:0
        for x = 0:(Nx-1)
            i = x+y*Nx;
            spin = bitget(s,i+1,numtype);
            fprintf('%d ',spin);
        end 
        fprintf('\n')
    end
    fprintf('\n')
end