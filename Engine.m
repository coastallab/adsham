clear all; close all; clc
%% Read setup-files
filename = 'Convergence';
data = readtable(strcat('Examples\',filename,'\',filename,'.txt'));
[epsilon, TWO_THETA, Corr, corrN, dt, width, length, nx, ny, dx, dy, T, nt, g, D, Inifile] = setParams(table2array(data(:,2)));
Bbar = importdata(strcat('Examples\',filename,'\','bathy.txt'));
[B, BX, BY] = reconstB(Bbar,nx,ny,width,length);

%% CFL check
CFL = sqrt(-g*min(min(Bbar))) * dt / min(dx,dy)

%% Boundary type
NorthBoundary = 'wall';
EastBoundary = 'periodic';
SouthBoundary = 'wall';
WestBoundary = 'periodic';

%% Variables
x = 0:dx:width; y=0:dy:length; [X,Y] = meshgrid(x,y);
U = zeros(ny+4,nx+4,4); % average value of w, hu, hv, hc
U0 = zeros(ny+4,nx+4,4); % average value of w, hu, hv, hc
hu = zeros(ny+3,nx+3,4); % reconstructed value of hu
hv = zeros(ny+3,nx+3,4); % reconstructed value of hv
hc = zeros(ny+3,nx+3,4); % reconstructed value of hc
w = zeros(ny+4,nx+4,4); % reconstructed value of h
h = zeros(ny+4,nx+4,4); % reconstructed value of h
u = zeros(ny+4,nx+4,4); % reconstructed value of u
v = zeros(ny+4,nx+4,4); % reconstructed value of v
c = zeros(ny+4,nx+4,4); % reconstructed value of c
dbydt = zeros(ny+4,nx+4,4,4); % cell value of temporal change of conservative variables
xflux = zeros(ny+4,nx+4,4); % reconstructed numerical flux along the x-direction
yflux = zeros(ny+4,nx+4,4); % reconstructed numerical flux along the x-direction
aplus = zeros(ny+2,nx+2); aminus = zeros(ny+2,nx+2); bplus = zeros(ny+2,nx+2); bminus = zeros(ny+2,nx+2);

%% Initial conditions
if Inifile == 1
    load(strcat('Examples\',filename,'\','Ini.mat'));
    IniW = Ini(:,:,1); IniHU = Ini(:,:,2); IniHV = Ini(:,:,3); IniHC = Ini(:,:,4);
    U0(3:ny+2,3:nx+2,1) = IniW;
    U0(3:ny+2,3:nx+2,2) = IniHU;
    U0(3:ny+2,3:nx+2,3) = IniHV;
    U0(3:ny+2,3:nx+2,4) = IniHC;
end

%% Numerical computation
count = 0;
U0 = BoundaryCondition(nx,ny,U0,NorthBoundary,EastBoundary,SouthBoundary,WestBoundary,0);
for k = 1:round(T/dt)
    if k < 3
        timeScheme = 1; % Euler
    elseif Corr == 0
        timeScheme = 2; % Adams-Bashforth
    elseif Corr == 1
        timeScheme = 3; % Adams-Moulton
    end
    
    [w(2:ny+3,2:nx+3,4),w(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,1),U0(2:ny+3,2:nx+3,1),U0(2:ny+3,3:nx+4,1),TWO_THETA);
    [w(2:ny+3,2:nx+3,3),w(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,1),U0(2:ny+3,2:nx+3,1),U0(3:ny+4,2:nx+3,1),TWO_THETA);
    [hu(2:ny+3,2:nx+3,4),hu(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,2),U0(2:ny+3,2:nx+3,2),U0(2:ny+3,3:nx+4,2),TWO_THETA);
    [hu(2:ny+3,2:nx+3,3),hu(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,2),U0(2:ny+3,2:nx+3,2),U0(3:ny+4,2:nx+3,2),TWO_THETA);
    [hv(2:ny+3,2:nx+3,4),hv(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,3),U0(2:ny+3,2:nx+3,3),U0(2:ny+3,3:nx+4,3),TWO_THETA);
    [hv(2:ny+3,2:nx+3,3),hv(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,3),U0(2:ny+3,2:nx+3,3),U0(3:ny+4,2:nx+3,3),TWO_THETA);
    [hc(2:ny+3,2:nx+3,4),hc(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,4),U0(2:ny+3,2:nx+3,4),U0(2:ny+3,3:nx+4,4),TWO_THETA);
    [hc(2:ny+3,2:nx+3,3),hc(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,4),U0(2:ny+3,2:nx+3,4),U0(3:ny+4,2:nx+3,4),TWO_THETA);
    
    [w(2:ny+3,2:nx+3,4),w(2:ny+3,2:nx+3,2)] = correction(BX(2:ny+3,2:nx+3),BX(2:ny+3,3:nx+4),U0(2:ny+3,2:nx+3,1),w(2:ny+3,2:nx+3,4),w(2:ny+3,2:nx+3,2));
    [w(2:ny+3,2:nx+3,3),w(2:ny+3,2:nx+3,1)] = correction(BY(2:ny+3,2:nx+3),BY(3:ny+4,2:nx+3),U0(2:ny+3,2:nx+3,1),w(2:ny+3,2:nx+3,3),w(2:ny+3,2:nx+3,1));
    
    h(2:ny+3,2:nx+3,:) = max(w(2:ny+3,2:nx+3,:)-cat(3,BY(3:ny+4,2:nx+3),BX(2:ny+3,3:nx+4),BY(2:ny+3,2:nx+3),BX(2:ny+3,2:nx+3)),0);
    [u(2:ny+3,2:nx+3,:), v(2:ny+3,2:nx+3,:), c(2:ny+3,2:nx+3,:)] = CalcUVC(h(2:ny+3,2:nx+3,:), hu(2:ny+3,2:nx+3,:), hv(2:ny+3,2:nx+3,:), hc(2:ny+3,2:nx+3,:), epsilon);
    
    aplus(2:ny+2,2:nx+2)  = max(max(u(2:ny+2,2:nx+2,2) + sqrt(g*h(2:ny+2,2:nx+2,2)), u(2:ny+2,3:nx+3,4) + sqrt(g*h(2:ny+2,3:nx+3,4))),0);
    aminus(2:ny+2,2:nx+2) = min(min(u(2:ny+2,2:nx+2,2) - sqrt(g*h(2:ny+2,2:nx+2,2)), u(2:ny+2,3:nx+3,4) - sqrt(g*h(2:ny+2,3:nx+3,4))),0);
    bplus(2:ny+2,2:nx+2)  = max(max(v(2:ny+2,2:nx+2,1) + sqrt(g*h(2:ny+2,2:nx+2,1)), v(3:ny+3,2:nx+2,3) + sqrt(g*h(3:ny+3,2:nx+2,3))),0);
    bminus(2:ny+2,2:nx+2) = min(min(v(2:ny+2,2:nx+2,1) - sqrt(g*h(2:ny+2,2:nx+2,1)), v(3:ny+3,2:nx+2,3) - sqrt(g*h(3:ny+3,2:nx+2,3))),0);
    
    delta = 10^-12;
    phix(2:ny+2,2:nx+2) = CAntiDissipation(u(2:ny+2,3:nx+3,4), u(2:ny+2,2:nx+2,2), sqrt(g*h(2:ny+2,3:nx+3,4)), sqrt(g*h(2:ny+2,2:nx+2,2)), delta);
    phiy(2:ny+2,2:nx+2) = CAntiDissipation(v(3:ny+3,2:nx+2,3), v(2:ny+2,2:nx+2,1), sqrt(g*h(3:ny+3,2:nx+2,3)), sqrt(g*h(2:ny+2,2:nx+2,1)), delta);
    
    xflux(2:ny+2,2:nx+2,1) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
        h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4), ...
        h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2), ...
        h(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2));
    xflux(2:ny+2,2:nx+2,2) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
        h(2:ny+2,3:nx+3,4) .* (u(2:ny+2,3:nx+3,4).^2 + 0.5 * g * h(2:ny+2,3:nx+3,4)), ...
        h(2:ny+2,2:nx+2,2) .* (u(2:ny+2,2:nx+2,2).^2 + 0.5 * g * h(2:ny+2,2:nx+2,2)), ...
        h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2));
    xflux(2:ny+2,2:nx+2,3) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
        h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4) .* v(2:ny+2,3:nx+3,4), ...
        h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2) .* v(2:ny+2,2:nx+2,2), ...
        h(2:ny+2,3:nx+3,4) .* v(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2) .* v(2:ny+2,2:nx+2,2));
    xflux(2:ny+2,2:nx+2,4) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
        h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4) .* c(2:ny+2,3:nx+3,4), ...
        h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2) .* c(2:ny+2,2:nx+2,2), ...
        phix(2:ny+2,2:nx+2) .* (h(2:ny+2,3:nx+3,4) .* c(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2) .* c(2:ny+2,2:nx+2,2)));
    
    yflux(2:ny+2,2:nx+2,1) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
        h(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3), ...
        h(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1), ...
        h(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1));
    yflux(2:ny+2,2:nx+2,2) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
        h(3:ny+3,2:nx+2,3) .* u(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3), ...
        h(2:ny+2,2:nx+2,1) .* u(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1), ...
        h(3:ny+3,2:nx+2,3) .* u(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1) .* u(2:ny+2,2:nx+2,1));
    yflux(2:ny+2,2:nx+2,3) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
        h(3:ny+3,2:nx+2,3) .* (v(3:ny+3,2:nx+2,3).^2 + 0.5 * g * h(3:ny+3,2:nx+2,3)), ...
        h(2:ny+2,2:nx+2,1) .* (v(2:ny+2,2:nx+2,1).^2 + 0.5 * g * h(2:ny+2,2:nx+2,1)), ...
        h(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1));
    yflux(2:ny+2,2:nx+2,4) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
        h(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3) .* c(3:ny+3,2:nx+2,3), ...
        h(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1) .* c(2:ny+2,2:nx+2,1), ...
        phiy(2:ny+2,2:nx+2) .* (h(3:ny+3,2:nx+2,3) .* c(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1) .* c(2:ny+2,2:nx+2,1)));
    
    hbar(3:ny+2,3:nx+2) = max(U0(3:ny+2,3:nx+2,1)-B(3:ny+2,3:nx+2),0);
    
    hc_by_dx_dx(3:ny+2,3:nx+2) = D(1) * (U0(3:ny+2,4:nx+3,4) - 2 * U0(3:ny+2,3:nx+2,4) + U0(3:ny+2,2:nx+1,4)) / dx^2;
    hc_by_dx_dy(3:ny+2,3:nx+2) = 0.25 * D(2) * (U0(4:ny+3,4:nx+3,4) - U0(4:ny+3,2:nx+1,4) - U0(2:ny+1,4:nx+3,4) + U0(2:ny+1,2:nx+1,4)) / dx / dy;
    hc_by_dy_dy(3:ny+2,3:nx+2) = D(3) * (U0(4:ny+3,3:nx+2,4) - 2 * U0(3:ny+2,3:nx+2,4) + U0(2:ny+1,3:nx+2,4)) / dy^2;
    
    source_term(3:ny+2,3:nx+2,:) = cat(3, zeros(ny,nx), - g / dx * hbar(3:ny+2,3:nx+2) .* (BX(3:ny+2,4:nx+3) - BX(3:ny+2,3:nx+2)), - g / dy * hbar(3:ny+2,3:nx+2) .* (BY(4:ny+3,3:nx+2) - BY(3:ny+2,3:nx+2)), hc_by_dx_dx(3:ny+2,3:nx+2)+2*hc_by_dx_dy(3:ny+2,3:nx+2)+hc_by_dy_dy(3:ny+2,3:nx+2));
    dbydt(3:ny+2,3:nx+2,:,3) = squeeze(xflux(3:ny+2,2:nx+1,:)-xflux(3:ny+2,3:nx+2,:)) / dx + squeeze(yflux(2:ny+1,3:nx+2,:)-yflux(3:ny+2,3:nx+2,:)) / dy + source_term(3:ny+2,3:nx+2,:);
    
    if timeScheme == 1 % Euler
        U(3:ny+2,3:nx+2,:) = U0(3:ny+2,3:nx+2,:) + dt * dbydt(3:ny+2,3:nx+2,:,3);
    else % Predictor, Adams-Bashforth
        U(3:ny+2,3:nx+2,:) = U0(3:ny+2,3:nx+2,:) + dt / 12 * (23 * dbydt(3:ny+2,3:nx+2,:,3) - 16 * dbydt(3:ny+2,3:nx+2,:,1) + 5 * dbydt(3:ny+2,3:nx+2,:,2));
    end
    U = BoundaryCondition(nx,ny,U,NorthBoundary,EastBoundary,SouthBoundary,WestBoundary,k);
    
    if timeScheme == 3 % Corrector, Adams-Moulton
        Ucur = U;
        for m = 1:corrN
            U0 = U;
            
            [w(2:ny+3,2:nx+3,4),w(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,1),U0(2:ny+3,2:nx+3,1),U0(2:ny+3,3:nx+4,1),TWO_THETA);
            [w(2:ny+3,2:nx+3,3),w(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,1),U0(2:ny+3,2:nx+3,1),U0(3:ny+4,2:nx+3,1),TWO_THETA);
            [hu(2:ny+3,2:nx+3,4),hu(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,2),U0(2:ny+3,2:nx+3,2),U0(2:ny+3,3:nx+4,2),TWO_THETA);
            [hu(2:ny+3,2:nx+3,3),hu(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,2),U0(2:ny+3,2:nx+3,2),U0(3:ny+4,2:nx+3,2),TWO_THETA);
            [hv(2:ny+3,2:nx+3,4),hv(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,3),U0(2:ny+3,2:nx+3,3),U0(2:ny+3,3:nx+4,3),TWO_THETA);
            [hv(2:ny+3,2:nx+3,3),hv(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,3),U0(2:ny+3,2:nx+3,3),U0(3:ny+4,2:nx+3,3),TWO_THETA);
            [hc(2:ny+3,2:nx+3,4),hc(2:ny+3,2:nx+3,2)] = reconstruction(U0(2:ny+3,1:nx+2,4),U0(2:ny+3,2:nx+3,4),U0(2:ny+3,3:nx+4,4),TWO_THETA);
            [hc(2:ny+3,2:nx+3,3),hc(2:ny+3,2:nx+3,1)] = reconstruction(U0(1:ny+2,2:nx+3,4),U0(2:ny+3,2:nx+3,4),U0(3:ny+4,2:nx+3,4),TWO_THETA);
            
            [w(2:ny+3,2:nx+3,4),w(2:ny+3,2:nx+3,2)] = correction(BX(2:ny+3,2:nx+3),BX(2:ny+3,3:nx+4),U0(2:ny+3,2:nx+3,1),w(2:ny+3,2:nx+3,4),w(2:ny+3,2:nx+3,2));
            [w(2:ny+3,2:nx+3,3),w(2:ny+3,2:nx+3,1)] = correction(BY(2:ny+3,2:nx+3),BY(3:ny+4,2:nx+3),U0(2:ny+3,2:nx+3,1),w(2:ny+3,2:nx+3,3),w(2:ny+3,2:nx+3,1));
            
            h(2:ny+3,2:nx+3,:) = max(w(2:ny+3,2:nx+3,:)-cat(3,BY(3:ny+4,2:nx+3),BX(2:ny+3,3:nx+4),BY(2:ny+3,2:nx+3),BX(2:ny+3,2:nx+3)),0);
            [u(2:ny+3,2:nx+3,:), v(2:ny+3,2:nx+3,:), c(2:ny+3,2:nx+3,:)] = CalcUVC(h(2:ny+3,2:nx+3,:), hu(2:ny+3,2:nx+3,:), hv(2:ny+3,2:nx+3,:), hc(2:ny+3,2:nx+3,:), epsilon);
            
            aplus(2:ny+2,2:nx+2)  = max(max(u(2:ny+2,2:nx+2,2) + sqrt(g*h(2:ny+2,2:nx+2,2)), u(2:ny+2,3:nx+3,4) + sqrt(g*h(2:ny+2,3:nx+3,4))),0);
            aminus(2:ny+2,2:nx+2) = min(min(u(2:ny+2,2:nx+2,2) - sqrt(g*h(2:ny+2,2:nx+2,2)), u(2:ny+2,3:nx+3,4) - sqrt(g*h(2:ny+2,3:nx+3,4))),0);
            bplus(2:ny+2,2:nx+2)  = max(max(v(2:ny+2,2:nx+2,1) + sqrt(g*h(2:ny+2,2:nx+2,1)), v(3:ny+3,2:nx+2,3) + sqrt(g*h(3:ny+3,2:nx+2,3))),0);
            bminus(2:ny+2,2:nx+2) = min(min(v(2:ny+2,2:nx+2,1) - sqrt(g*h(2:ny+2,2:nx+2,1)), v(3:ny+3,2:nx+2,3) - sqrt(g*h(3:ny+3,2:nx+2,3))),0);
            
            delta = 10^-12;
            phix(2:ny+2,2:nx+2) = CAntiDissipation(u(2:ny+2,3:nx+3,4), u(2:ny+2,2:nx+2,2), sqrt(g*h(2:ny+2,3:nx+3,4)), sqrt(g*h(2:ny+2,2:nx+2,2)), delta);
            phiy(2:ny+2,2:nx+2) = CAntiDissipation(v(3:ny+3,2:nx+2,3), v(2:ny+2,2:nx+2,1), sqrt(g*h(3:ny+3,2:nx+2,3)), sqrt(g*h(2:ny+2,2:nx+2,1)), delta);
            
            xflux(2:ny+2,2:nx+2,1) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
                h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4), ...
                h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2), ...
                h(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2));
            xflux(2:ny+2,2:nx+2,2) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
                h(2:ny+2,3:nx+3,4) .* (u(2:ny+2,3:nx+3,4).^2 + 0.5 * g * h(2:ny+2,3:nx+3,4)), ...
                h(2:ny+2,2:nx+2,2) .* (u(2:ny+2,2:nx+2,2).^2 + 0.5 * g * h(2:ny+2,2:nx+2,2)), ...
                h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2));
            xflux(2:ny+2,2:nx+2,3) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
                h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4) .* v(2:ny+2,3:nx+3,4), ...
                h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2) .* v(2:ny+2,2:nx+2,2), ...
                h(2:ny+2,3:nx+3,4) .* v(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2) .* v(2:ny+2,2:nx+2,2));
            xflux(2:ny+2,2:nx+2,4) = HLLsolver(aplus(2:ny+2,2:nx+2), aminus(2:ny+2,2:nx+2), ...
                h(2:ny+2,3:nx+3,4) .* u(2:ny+2,3:nx+3,4) .* c(2:ny+2,3:nx+3,4), ...
                h(2:ny+2,2:nx+2,2) .* u(2:ny+2,2:nx+2,2) .* c(2:ny+2,2:nx+2,2), ...
                h(2:ny+2,3:nx+3,4) .* c(2:ny+2,3:nx+3,4) - h(2:ny+2,2:nx+2,2) .* c(2:ny+2,2:nx+2,2));
            
            yflux(2:ny+2,2:nx+2,1) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
                h(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3), ...
                h(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1), ...
                h(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1));
            yflux(2:ny+2,2:nx+2,2) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
                h(3:ny+3,2:nx+2,3) .* u(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3), ...
                h(2:ny+2,2:nx+2,1) .* u(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1), ...
                h(3:ny+3,2:nx+2,3) .* u(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1) .* u(2:ny+2,2:nx+2,1));
            yflux(2:ny+2,2:nx+2,3) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
                h(3:ny+3,2:nx+2,3) .* (v(3:ny+3,2:nx+2,3).^2 + 0.5 * g * h(3:ny+3,2:nx+2,3)), ...
                h(2:ny+2,2:nx+2,1) .* (v(2:ny+2,2:nx+2,1).^2 + 0.5 * g * h(2:ny+2,2:nx+2,1)), ...
                h(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1));
            yflux(2:ny+2,2:nx+2,4) = HLLsolver(bplus(2:ny+2,2:nx+2), bminus(2:ny+2,2:nx+2), ...
                h(3:ny+3,2:nx+2,3) .* v(3:ny+3,2:nx+2,3) .* c(3:ny+3,2:nx+2,3), ...
                h(2:ny+2,2:nx+2,1) .* v(2:ny+2,2:nx+2,1) .* c(2:ny+2,2:nx+2,1), ...
                h(3:ny+3,2:nx+2,3) .* c(3:ny+3,2:nx+2,3) - h(2:ny+2,2:nx+2,1) .* c(2:ny+2,2:nx+2,1));
            
            hbar(3:ny+2,3:nx+2) = max(U0(3:ny+2,3:nx+2,1)-B(3:ny+2,3:nx+2),0);
            source_term(3:ny+2,3:nx+2,:) = cat(3, zeros(ny,nx), - g / dx * hbar(3:ny+2,3:nx+2) .* (BX(3:ny+2,4:nx+3) - BX(3:ny+2,3:nx+2)), - g / dy * hbar(3:ny+2,3:nx+2) .* (BY(4:ny+3,3:nx+2) - BY(3:ny+2,3:nx+2)), zeros(ny,nx));
            dbydt(3:ny+2,3:nx+2,:,4) = squeeze(xflux(3:ny+2,2:nx+1,:)-xflux(3:ny+2,3:nx+2,:)) / dx + squeeze(yflux(2:ny+1,3:nx+2,:)-yflux(3:ny+2,3:nx+2,:)) / dy + source_term(3:ny+2,3:nx+2,:);
            U(3:ny+2,3:nx+2,:) = Ucur(3:ny+2,3:nx+2,:) + dt / 24 * (9 * dbydt(3:ny+2,3:nx+2,:,4) + 19 * dbydt(3:ny+2,3:nx+2,:,3) - 5 * dbydt(3:ny+2,3:nx+2,:,1) + 1 * dbydt(3:ny+2,3:nx+2,:,2));
            U = BoundaryCondition(nx,ny,U,NorthBoundary,EastBoundary,SouthBoundary,WestBoundary,k);
        end
    end
    U0 = U; % current state to old state
    dbydt = dbydt(:,:,:,[3,1,2]); % [t-2, t-3, t-1] -> [t-1, t-2, t-3]
    
    % Visualization
    if mod(k,50) == 0
        plot(x,U(4,3:nx+2,1))
        title(['Time = ',num2str(k*dt),' s'])
        shading flat
        pause(0.001)
    end
    if mod(k,T/dt)==0
        count = count + 1;
        output = squeeze(U(3:ny+2,3:nx+2,:));
    end
end
