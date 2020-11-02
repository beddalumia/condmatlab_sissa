function Ts = cyclebits(s,N,numtype)
    oflowbit = bitget(s,N,numtype);
    cutted = bitset(s,N,0,numtype);
    shifted = bitshift(cutted,1,numtype);
    Ts = bitset(shifted,1,oflowbit,numtype);
end