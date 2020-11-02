N=16
load('LevelsChain16.mat') % Loads: E = cell(1,16) :)
AllData = [];
for i = 1:N
y = real(E{i});
AllData = [AllData, y'];
end

hfit = histfit(AllData-min(AllData),2*N,'poisson')
%h = histogram(AllData./N,2*N,'Normalization','pdf');