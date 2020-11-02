function counter = countSzKstates(vec,k,N,numtype) 
    M = length(vec);
    counter = 0;
    for a = 1:M
        R = checkstate(vec(a),k,N,numtype);
        if R >= 0
            counter = counter + 1;
        end
    end
   end