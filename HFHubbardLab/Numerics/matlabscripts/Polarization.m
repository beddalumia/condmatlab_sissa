clear all
a = linspace(-1,1,250);
t = linspace(-1,1,250);

for i = 1:length(a)
    for j = 1:length(t)
        x = a(i);
        y = t(j);
        P(i,j) = (2*x*(sqrt(x^2+y^2)-y))/(x^2+(sqrt(x^2+y^2)-y)^2);
    end
end
s = surf(a,t,P); hold on
s.EdgeColor = 'none';
s.FaceAlpha = 0.5;
plot3(zeros(1,length(a)),a,sign(a))

%print(gcf,'PolarizationTerm','-dpng','-r1920')