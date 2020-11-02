N=16;
load('LevelsChain16.mat') % Loads: E = cell(1,16) :)

%% We need to sort all eigenvalues 

AllData = [];
for i = 1:N
y = real(E{i});
AllData = [AllData, y'];
end

Levels = sort(AllData);

%% Let's define the Spacing-list 

Spacing = diff(Levels);
Spacing = Spacing(Spacing>3*10^(-14)); % Here the cutoff

%% And so the Phys. Rev. B 75, 155111 (2007) Statistics 

M = length(Spacing)-1;
Ratio = -ones(1,M);
for j = 1:M
    A = Spacing(j);
    B = Spacing(j+1);
    Ratio(j) = min(A,B) / max(A,B);
end

h = histogram(Ratio,20,'Normalization','pdf');
y = h.Values;
x = h.BinEdges;
delta = x(2)-x(1);
x = x+delta;
x = x(1:length(x)-1);
%plot(x,y,'.-');
xx = linspace(0,1,100);
theory = 2./(1+xx).^2;
hold on
plot(xx,theory,'-r');
xlim([0,1])