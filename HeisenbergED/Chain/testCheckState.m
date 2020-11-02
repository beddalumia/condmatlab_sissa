N = 8;
nbit = 32;
numtype = sprintf('int%d',nbit);

clc
for k = -N/2+1:N/2
    counter = 0;
    fprintf('k = %f, \n', k);
    for s = 0:2^N-1
        R = checkstate(s,k,N,numtype);
        if R >= 0
            %fprintf('\nState: |%s> (= %d), ',dec2bin(s,8), s);
            %fprintf('R = %d \n', R);
            counter = counter + 1;
        end
    end
    fprintf('counter_k = %d \n\n', counter);
end