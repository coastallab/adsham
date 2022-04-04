function [u, v, c] = CalcUVC(h, hu, hv, hc, epsilon)
divide_by_h = sqrt(2) * h ./ sqrt(h.^4 + max(h.^4, epsilon));
u = divide_by_h .* hu;
v = divide_by_h .* hv;
c = divide_by_h .* hc;
end