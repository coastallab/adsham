function phi = CAntiDissipation(uL, uR, aL, aR, delta)
Fr = (0.5 * (1 + sign(abs(uL)-abs(uR))) .* abs(uL) ./ aL + 0.5 * (1 - sign(abs(uL)-abs(uR))) .* abs(uR) ./ aR) .* and(aL~=0,aR~=0);
phi = (Fr + delta) ./ (Fr + 1);
end
