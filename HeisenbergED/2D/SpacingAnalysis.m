Nx = 4;
Ny = 4;
N = Nx*Ny;
Name = sprintf('Levels%dx%d.mat',Nx,Ny);
load(Name); % Loads: E ~ cell(Nx,Ny) :)
Levels = E;
AllData = [];
for i = 1%:Nx
 for j = 1%:Ny
    y = real(Levels{i,j})/N;
    AllData = [AllData; y];
 end
end

%% We need to sort all eigenvalues 
Levels = sort(AllData);

%% Let's define the Spacing-list 
Spacing = diff(Levels);
Spacing = Spacing(Spacing>3*10^(-5)); % Here the cutoff
%Spacing = Spacing(Spacing~=0);
%% And so the Phys. Rev. B 75, 155111 (2007) Statistics 

M = length(Spacing)-1;
Ratio = -ones(1,M);
for j = 1:M
    A = Spacing(j);
    B = Spacing(j+1);
    Ratio(j) = min(A,B) / max(A,B);
end

xx = linspace(0,1,100);
theory = 2./(1+xx).^2;

figure();
h = histogram(Ratio,20,'Normalization','pdf'); hold on
plot(xx,theory,'-r','Linewidth',2); hold off
legend('Binned Data','Th. Prediction')
title(sprintf('%d x %d', Nx, Ny))
xlim([0,1])

% figure("Name",'curve');
% y = h.Values;
% x = h.BinEdges;
% delta = x(2)-x(1);
% x = x+delta;
% x = x(1:length(x)-1);
% plot(x,y); hold on
% plot(xx,theory,'-r');
% xlim([0,1])
