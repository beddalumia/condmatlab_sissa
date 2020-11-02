%% Importing Data... 
Nx = 2;
Ny = 10;
N = Nx*Ny;
Name = sprintf('AllMagnetizations%dx%d.mat',Nx,Ny);
load(Name); % Loads: E ~ cell(1,N+1) :)
Levels = E;
AllData = [];
gs = 0;
for i = 1:N+1
    y = real(Levels{i});
    if gs >= min(y)
        gs = min(y);
        IGS = i;
    end
    AllData = [AllData; y];
end

%% Density of States
% figure("Name",'DOS')
% DOS = histogram(AllData./N,20,'Normalization','pdf');
figure("Name",'DOS')
DOSfit = histfit(AllData-min(AllData),20,'poisson');


%% We need to sort all eigenvalues 
Levels = sort(AllData)./N;

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

figure("Name",'Spacing Analysis')
h = histogram(Ratio,20,'Normalization','pdf');
y = h.Values;
x = h.BinEdges;
delta = x(2)-x(1);
x = x+delta;
x = x(1:length(x)-1);
xx = linspace(0,1,100);
theory = 2./(1+xx).^2;
hold on
plot(xx,theory,'-r','LineWidth',2);
xlim([0,1])