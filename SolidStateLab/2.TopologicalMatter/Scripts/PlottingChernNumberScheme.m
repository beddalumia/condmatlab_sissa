N = 7;
figure("Name",'HH_Lattice')

for i=1:N
    xline(i,'k'); hold on
    yline(i,'k');
    for j = 1:N
        if i < N && j < N
            plot(i,j,'o','MarkerSize',6,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
%             plot(i+0.5,j+0.5,'.','MarkerSize',6,'MarkerEdgeColor','r');
%             plot(i+0.5,j+0.5,'o','MarkerSize',6,'MarkerEdgeColor','r','LineWidth',1);
              text(i+0.40,j+0.55,'\gamma','Color',[0.9,i/10,j/10],'fontweight', 'bold')
        else
            plot(i,j,'o','MarkerSize',6,'MarkerEdgeColor','k', 'MarkerFaceColor', [0.7,0.7,0.7]);
        end
    end
end
axis equal
xlim([0.5,N+0.5]);
ylim([0.5,N+0.5]);
xticks([]);
yticks([]);
set(gca,'visible','off')