function [L,R] = reconstruction(UL, UC, UR, TWO_THETA)
dx_grad_over_two = 0.25 * minmod(TWO_THETA * (UC - UL), (UR - UL), TWO_THETA * (UR - UC));
R = UC + dx_grad_over_two;
L = UC - dx_grad_over_two;
end