function counter = countKstates(k,N,numtype)  
    counter = 0;
    for s = 0:2^N-1
        R = checkstate(s,k,N,numtype);
        if R >= 0
            counter = counter + 1;
        end
    end
   end