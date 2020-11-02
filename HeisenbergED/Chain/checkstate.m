function R = checkstate(s,k,N,numtype)
    R = -1; t = s;
    for i = 1:N
        t = cyclebits(t,N,numtype);
        if t < s % It's sure that s is not the representative
           return;
        elseif t == s % Round trip ended: s could be the representative
            if mod(k,N/i) == 0 % Compatible with k (eq. 119)
            R = i; return;
            else
                return; % Not compatible
            end
        end
    end
end

% There is an implicit "if t > s the roundtrip is still going so we don't
% know anything..."