Nx = 1;
Ny = 20;
N = Nx*Ny;
Name = sprintf('Levels%dx%d.mat',Nx,Ny);
load(Name); % Loads: E ~ cell(Nx,Ny) :)
Levels = E;
FullBloch = [];
for i = 1:Nx
 for j = 1:Ny
    y = real(Levels{i,j})/N;
    FullBloch = [FullBloch;y];
 end
end
FullBloch = sort(FullBloch);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N == 16                                      %
Check = sprintf('Check%dx%d.mat',Nx,Ny);        %
load(Check) % Mz=0 Diag for check ;)            %
ZeroMagn = E/N;                                 %
end                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Nx ~= Ny && Nx ~= 1 && N == 16
    emaN = sprintf('Levels%dx%d.mat',Ny,Nx);
    load(emaN); % Loads: E ~ cell(Ny,Nx) :)
    sleveL = E;

BlochFull = [];
for i = 1:Ny
 for j = 1:Nx
    y = real(sleveL{i,j})/N;
    BlochFull = [BlochFull; y];
 end
end
BlochFull = sort(BlochFull);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nx ~= Ny && Nx~=1 && N == 16                       %
    checK = sprintf('Check%dx%d.mat',Ny,Nx);          %
    load(checK) % Mz=0 Diag for check ;)              %
    MagnZero = E/N;                                   %
end                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure("Name",'DOSfit')
hfit = histfit(FullBloch*N-min(FullBloch*N),2*N,'poisson');
legend('Binned Data', 'Fitted Poisson','Location','best')
figure("Name", 'DOS')
h1 = histogram(FullBloch,50,'Normalization','pdf'); hold on
bins = h1.BinEdges;
dos = smooth(h1.Values);
bins = bins(2:length(bins));
delta = bins(2)-bins(1);
energies = bins-delta/2;
plot(energies,dos,'Linewidth',2);
legend('Binned Data', 'Smoothed Bins','Location','best')
%h2 = histogram(ZeroMagn,50,'Normalization','pdf','BinEdges',bins);

legend();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARISONS                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N == 16
flag = 0;
for i = 1:length(FullBloch)
    if FullBloch(i) - ZeroMagn(i) > 10^(-12)
        flag = flag+1;
    end
end
if flag > 0
    fprintf('ERROR: %d eigenvalues are different\n', flag);
end

figure("Name",'Comparisons')
plot(1:length(FullBloch),FullBloch,'-','LineWidth',2)
hold on
plot(1:length(FullBloch),ZeroMagn,'--','LineWidth',2)
if Nx ~= Ny && Nx~=1
    plot(1:length(FullBloch),BlochFull,':','LineWidth',1);
    plot(1:length(FullBloch),MagnZero,':','LineWidth',1);
end
end