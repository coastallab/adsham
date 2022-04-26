%% Result display
figure
pcolor(X,Y,result(:,:,4,1))
shading flat
xlabel('x (m)'), ylabel('y (m)')
c = colorbar
c.Label.String = 'C';
title('Time = 7.5 s')
dir = ['../../../results/D',num2str(D(1)),'T7.5s','.png'];
saveas(gcf,dir)
close 

figure
pcolor(X,Y,result(:,:,4,1))
shading flat
xlabel('x (m)'), ylabel('y (m)')
c = colorbar
c.Label.String = 'C';
title('Time = 15 s')
dir = ['../../../results/D',num2str(D(1)),'T15s','.png'];
saveas(gcf,dir)
close 