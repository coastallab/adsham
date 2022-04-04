function [B, BX, BY] = reconstB(Bbar,nx,ny,width,length)
% B = zeros(ny,nx);
% BX = zeros(ny+1,nx+1);
% BY = zeros(ny+1,nx+1);
[x,y] = meshgrid(linspace(0,width,size(Bbar,2)),linspace(0,length,size(Bbar,1)));
[x1,y1] = meshgrid(linspace(-2.5*width/(nx-1),width+2.5*width/(nx-1),nx+5),linspace(-2.5*length/(ny-1),length+2.5*length/(ny-1),ny+5));
x1(x1<0) = 0; y1(y1<0) = 0; x1(x1>width) = width; y1(y1>length) = length;
B0 = interp2(x,y,Bbar,x1,y1);
for i = 1:ny+4
    for j = 1:nx+4
        B(i,j) = 0.25 * (B0(i,j) + B0(i+1,j) + B0(i,j+1) + B0(i+1,j+1));
        BX(i,j) = 0.5 * (B0(i+1,j) + B0(i,j));
        BY(i,j) = 0.5 * (B0(i,j) + B0(i,j+1));
    end
end
BX(:,end+1) = 0.5 * (B0(2:end,end) + B0(1:end-1,end));
BY(end+1,:) = 0.5 * (B0(end,2:end) + B0(end,1:end-1));
end