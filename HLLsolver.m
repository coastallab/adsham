function FLUX = HLLsolver(aplus, aminus, Fplus, Fminus, Udiff)
FLUX = (aplus .* Fminus - aminus .* Fplus + aplus .* aminus .* Udiff) ./ (aplus - aminus) .* (aplus - aminus ~= 0);
end