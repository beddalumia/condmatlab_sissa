h = findobj(gca,'Type','line')
x = h.XData;
y = h.YData;
Emin = min(y)
OptimalLatticeParameter = x(y==min(y))
hold on
scatter(OptimalLatticeParameter,Emin,125,'x','r','LineWidth',2);
xline(OptimalLatticeParameter,'--r');
yline(Emin,'--r');
axis square
xlabel('Lattice Parameter [Å]','Interpreter','tex');
ylabel('Cohesion Energy [cm^{-1}]','Interpreter','tex');
legend({'Computed Data', 'Polynomial Fit', sprintf('(%f,%f)',OptimalLatticeParameter,Emin)},'FontSize', 14) 
title('Geometry Optimization for HCP lattice')
set(gca,'FontSize',14)