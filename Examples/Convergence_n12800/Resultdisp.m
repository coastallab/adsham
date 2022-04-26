%% Result display
figure;
plot(linspace(0,1,length(x)),result(:,1))
xlabel('x (m)')
ylabel('water level (m)')
legend('Reference','location','northwest','box','off')
dir = ['../../../results/reference.png'];
saveas(gcf,dir)

