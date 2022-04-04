function [epsilon, TWO_THETA, Corr, corrN, dt, width, length, nx, ny, dx, dy, T, nt, g, D, Inifile] = setParams(data)
epsilon = data(1)^4;
TWO_THETA = data(2) * 2;
Corr = data(3);
corrN = data(4);
dt = data(5);
width = data(6);
length = data(7);
nx = data(8);
dx = width / (nx-1);
ny = data(9);
dy = length / (ny-1);
T = data(10);
nt = T / dt;
g = data(11);
D = [data(12), data(13), data(14)];
Inifile = data(15);
end
