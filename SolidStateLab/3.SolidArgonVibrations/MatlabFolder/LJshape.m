% Equilibrium Distance for diatomic molecule (experimental)
r0 = 3.758; %Å
% => sigma
sigma = r0/2^(1/6);
% Binding Energy (experimental)
e = 99.55; %cm^{-1}
% Some coordinates...
r = linspace(3,9,1000);
% LJ computation
U = LJpot(r);
%Plot
plot(r,U,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5); hold on
xline(sigma,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1); 
yline(-e,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1);
xline(r0,':')
yline(0,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1);
ylim([-110,110])
title('Lennard Jones model for Argon')
xlabel('Interatomic Distance [Å]')
ylabel('Interaction Energy [cm^{-1}]','interpreter','tex')
legend('$V_{\mathrm{LJ}}$','$\enspace\sigma$','$\enspace\varepsilon$','$\enspace r_0$','interpreter','latex')