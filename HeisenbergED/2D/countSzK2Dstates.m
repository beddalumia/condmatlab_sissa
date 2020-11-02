function counter = countSzK2Dstates(vec,kx,ky,Nx,Ny,numtype) 
    M = length(vec);
    counter = 0;
    for a = 1:M
        Ns = check2Dstate(vec(a),kx,ky,Nx,Ny,numtype);
        if Ns > 0
            counter = counter + 1;
        end
    end
   end