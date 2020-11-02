function Ts = Ty(s,Nx,Ny,numtype)
    oflowbit = zeros(1,Nx);
    cutted = s;
    for x = 1:Nx
        oflowbit(x) = bitget(s,x+(Ny-1)*Nx,numtype);
        cutted = bitset(cutted,x+(Ny-1)*Nx,0,numtype);
    end
    shifted = bitshift(cutted,Nx,numtype);
    Ts = shifted;
    for x = 1:Nx
        Ts = bitset(Ts,x,oflowbit(x),numtype);
    end
end