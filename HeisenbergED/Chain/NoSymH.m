N = 8;
H = zeros(2^N,2^N);
for a = int32(0:2^N-1)
    for i = 0:(N-1)
        j = mod(i+1,N);
        ai = bitget(a,i+1,'int32');
        aj = bitget(a,j+1,'int32');
        if ai == aj
            H(a+1,a+1) = H(a+1,a+1) + 1/4;
        else
            H(a+1,a+1) = H(a+1,a+1) - 1/4;
            b = bitxor(a,2^i+2^j,'int32');
            H(a+1,b+1) = 1/2;
        end
    end   
end
E = eig(H);