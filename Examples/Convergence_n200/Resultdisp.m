%% Result display
load('reference.mat')
ratio = 12800 / (length(x)-1);
werr = sum(abs(result(:,1)-reference(1:ratio:end,1)))/sum(abs(reference(1:ratio:end,1)));
huerr = sum(abs(result(:,2)-reference(1:ratio:end,2)))/sum(abs(reference(1:ratio:end,2)));
hcerr = sum(abs(result(:,4)-reference(1:ratio:end,4)))/sum(abs(reference(1:ratio:end,4)));

fprintf('Table. L1 error of w, hu and hc for given grid number N \n')
fprintf('||   N  ||    w     ||    hu    ||    hc    || \n')
fprintf('|| %4d || %.2e || %.2e || %.2e || \n',length(x)-1, werr, huerr, hcerr)

figure;
plot(linspace(0,1,length(x)),result(:,1))
hold on
plot(linspace(0,1,size(reference,1)),reference(:,1),'k--')
xlabel('x (m)')
ylabel('water level (m)')
str = ['N=',num2str(length(x)-1)];
legend(str,'Reference','location','northwest','box','off')
dir = ['../../../results/Convergence_N',num2str(length(x)-1),'.png'];
saveas(gcf,dir)

