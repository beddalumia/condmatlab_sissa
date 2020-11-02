% Physical System
t = 1;
ell = 1;
N = 2*10^2;

% Staggered Potential Values
a = linspace(0,3*t,100);

% Folded Brillouin Zone
dk = 2*pi/N;
k = linspace(dk,pi/ell,N/2);

% Imperturbed Band Structure
ek = -2*t*cos(k.*ell);

% Perturbed Band-Structure
Ek = zeros(length(a),N/2);
for i = 1:length(a)
        Ek(i,:) = sqrt(ek.^2 + (a(i))^2);
end

% Displaying
for j = 1:N/2
    plot(a,-Ek(:,j),'Color',[155/255, 194/255, 177/255]); hold on
    plot(a,+Ek(:,j),'Color',[240/255,230/255,140/255]); hold on
end
xlabel('Staggered Potential Strength [$a$]','Interpreter','latex')
ylabel('Single Particle Energies [$\varepsilon_k$]','Interpreter','latex')
title('Noninteracting Band-Structure [$N=200,\,t=1$]','Interpreter','latex')