function [r,L] = representative(s,N,numtype)
    r = s;  %
    t = s;  % Inizializations
    L = 0;  %
    for i = 1:N-1
        t = cyclebits(t,N,numtype);
        if t < r
           r = t;
           L = i;
        end
    end
end