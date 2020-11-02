a = 2.8; % Lattice parameter
b_1 = [1,0].*a;              %%%%%%%%%%%%%%%%%%%%%%
b_2 = [-0.5, sqrt(3)/2]*a;   % Bond displacements %
b_3 = [-0.5, - sqrt(3)/2]*a; %%%%%%%%%%%%%%%%%%%%%%

eps = 0;%-10.7 eV to vacuum but since eps = EF and we want to set EF = 0...
t = 2.5; %eV So we'll keep the energies in electronvolt

kmin = -pi/a; % Times 10 to calculate DOS
kmax = pi/a;  % Times 10 to calculate DOS
step = 0.01;
number_of_ks= ceil((kmax - kmin)/step);
[k_1,k_2] = meshgrid(kmin:step:kmax,kmin:step:kmax);

Scos = cos(k_1*b_1(1)+k_2*b_1(2))+cos(k_1*b_2(1)+ k_2*b_2(2))+cos(k_1*b_3(1)+k_2*b_3(2));
Ssin = sin(k_1*b_1(1)+k_2*b_1(2))+sin(k_1*b_2(1)+ k_2*b_2(2))+sin(k_1*b_3(1)+k_2*b_3(2));

Ec = eps + t*sqrt(Scos.^2 + Ssin.^2);
Ev = eps - t*sqrt(Scos.^2 + Ssin.^2);

figure('name','Band Structure of Graphene')
surf(k_1,k_2,Ec);
hold on
surf(k_1,k_2,Ev)
colormap summer
hold off
%print(gcf,'Band Structure','-bestfit','-dpdf','-r720')
%print(gcf,'Band Structure','-bestfit','-dpdf','painters')

%%%%%%%%%%%%%%%%%%%%%
% Density of States %
%%%%%%%%%%%%%%%%%%%%%
number_of_bins = 50;

% Counting the levels in a energy bin (E, E+dE)
[N_v,edgesv] = histcounts(Ev(:),number_of_bins);
[N_c,edgesc] = histcounts(Ec(:),number_of_bins);

% Normalization
D_v = N_v * (number_of_ks)^(-2);
D_c = N_c * (number_of_ks)^(-2);
N_tot = (sum(D_c) + sum(D_v)) % Check: has to be equal to 2 (#Bands)

% Bin centers
C = zeros(2,number_of_bins);
for i = 1:number_of_bins 
    C(1,i) = (edgesv(i)+edgesv(i+1))/2;
    C(2,i) = (edgesc(i)+edgesc(i+1))/2;
end

figure('name','Density of States of Graphene')
plot(C(1,1:number_of_bins),D_v, '.')
hold on
plot(C(2,1:number_of_bins),D_c, '.')
%print(gcf,'DOS','-bestfit','-dpdf','-painters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Smearing of the DOS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp = 0.1;     % Smearing Parameter

E = linspace(-8,8,500);
g_v = SMEARING(D_v,sp,C(1,:),E);
g_c = SMEARING(D_c,sp,C(2,:),E);

plot(E,g_v+g_c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized Gaussian 
function y = Gauss(x,mean,dev)
    y = 1/(sqrt(2*pi)*dev)*exp(-1/2*((x-mean)/dev).^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Smearing
function y = SMEARING(v,dev,centers,domain)
    y = zeros(1,length(domain));
    for i = 1:length(domain)
       for j = 1:length(v)
            y(i)= y(i) + v(j)*Gauss(domain(i),centers(j),dev)/length(v);
       end  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%