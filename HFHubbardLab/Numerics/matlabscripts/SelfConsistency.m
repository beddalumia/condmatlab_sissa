%% PHYSICAL SETUP %%
clear all
clc
% Imperturbed Band-Structure
t = 1;
ell = 1;
N = 10^4;
dk = 2*pi/N;
k = linspace(dk,pi/ell,N/2);
tk = -2*t*cos(k.*ell);

% Perturbative Physical Scales
  U = 4;%linspace(0.00001,1,250);
  a = linspace(-4,4,250);

%% MAIN %%

% Cycle
D = zeros(length(U),length(a));
M = zeros(length(U),length(a));
Emb = zeros(length(U),length(a));
Ek = zeros(length(U),length(a),N/2);
Eku = zeros(length(U),length(a),N/2);
Ekd = zeros(length(U),length(a),N/2);
for i = 1:length(U)
    for j = 1:length(a)
        
        % INIT VALUES for mean-field parameters
        Delta = 0.5;
        Magn = 0.5; 
        x = U(i)*(Delta+Magn)/2;
        y = U(i)*(Delta-Magn)/2;
        
        % Looking for Self-Consistency
        [X,Y] = selfconsistence(x,y,U(i),a(j),N,tk);
        D(i,j) = (X+Y)/U(i);
        M(i,j) = (X-Y)/U(i);
        Emb(i,j) = -X*Y/U(i);
        Eku(i,j,:) = sqrt(tk.^2 + (a(j)-X)^2);
        Ekd(i,j,:) = sqrt(tk.^2 + (a(j)-Y)^2);
        Ek(i,j,:) = Eku(i,j,:)+Ekd(i,j,:);
    end
end

% figure("Name",sprintf('Mean-Field Parameters Vs U [a fixed to %f]',a(10)));
% plot(U,D(:,10),'red'); hold on
% plot(U,M(:,10),'blue');
% 
% figure("Name",sprintf('Mean-Field Parameters Vs a [U fixed to %f]',U(1)));
% plot(a,D(1,:),'red'); hold on
% plot(a,M(1,:),'blue');
% 
% figure("Name",'Whole Parameter Space')
% surf(a,U,D); hold on
% surf(a,U,M);
% xlabel('a');
% ylabel('U');

figure("Name",'Magnetization and Staggering')
plot(a,M,'.-','Color',[0, 0.4470, 0.7410],'MarkerSize',6)
hold on
plot(a,D,'.-','Color',[0.8500, 0.3250, 0.0980],'MarkerSize',6)

figure("Name",'Mean-Field Single Particle Energies')
for ik = 1:ceil(N/50):N/2
        plot(a,Ek(1:length(U),1:length(a),ik),'Color',[240/255,230/255,140/255]); hold on
        plot(a,-Ek(1:length(U),1:length(a),ik),'Color',[155/255, 194/255, 177/255])
        plot(a,Emb(1:length(U),1:length(a)),'--r')
end

figure("Name",'Mean-Field Single Particle Energies [spin-polarized]')
for ik = 1:ceil(N/200):N/2
        plot(a,-Ekd(1:length(U),1:length(a),ik),'Color',[1 0 0 0.3]); hold on
        plot(a,-Eku(1:length(U),1:length(a),ik),'Color',[0 0 1 0.3]); 
        plot(a,Ekd(1:length(U),1:length(a),ik),'Color',[1 0 0 0.3]);
        plot(a,Eku(1:length(U),1:length(a),ik),'Color',[0 0 1 0.3]);
end

%% SELF CONSISTENCY ROUTINE %%

function [x,y] = selfconsistence(x,y,U,a,N,tk)

% Self-Consistency definition parameters
SELFx = 0.00001;
SELFy = 0.00001;

% Mixing Parameter -> Tunable Speed
SELFmix = 0.9;
PRODmix = 1-SELFmix;

stepCOUNTER = 0;
exitFLAG = false;
while exitFLAG ~= true

    % Computing new parameters
    X = U/N*sum((a-y)./sqrt(tk.^2 + (a-y).^2));
    Y = U/N*sum((a-x)./sqrt(tk.^2 + (a-x).^2));
    if imag(X) ~= 0 || imag(Y) ~= 0
       fprinft('ERROR')
       break
    end

    % Computing relative distances from self-consistency
    Dx = abs(x-X); dx = norm(Dx/x);
    Dy = abs(y-Y); dy = norm(Dy/y);
    
    % ~Comparison and Move~
    if (dx < SELFx && dy < SELFy) || stepCOUNTER > 10*N
        exitFLAG = true;
    else
        x = SELFmix * x + PRODmix * X;
        y = SELFmix * y + PRODmix * Y;
        %clc
%         fprintf('\nUPDATE of MEAN-FIELD Parameters!\n');
%         fprintf('********************************\n');
%         fprintf('Delta = %f\n', real(x+y)/U);
%         fprintf('Magn = %f\n', real(x-y)/U);
    end
    stepCOUNTER = stepCOUNTER + 1;
end
Delta = (x+y)/U;
Magn = (x-y)/U;
fprintf('\n********************************');
fprintf('\nConvergence achieved in %d steps\n', stepCOUNTER);
fprintf('********************************\n');
fprintf('Accuracy on Delta+Magn: %f\n', dx);
fprintf('Accuracy on Delta-Magn: %f\n', dy);
fprintf('********************************\n');
fprintf('********************************\n');
end