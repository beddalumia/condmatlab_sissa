function Ts = Tx(s,Nx,Ny,numtype)
    oflowbit = zeros(1,Ny);
    for y = 1:Ny
        oflowbit(y) = bitget(s,y*Nx,numtype);
    end
    cutted = bitset(s,Nx*Ny,0,numtype);
    shifted = bitshift(cutted,1,numtype);
    Ts = shifted;
    for y = 1:Ny
        Ts = bitset(Ts,(y-1)*Nx+1,oflowbit(y),numtype);
    end
end