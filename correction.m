function [wLcorr,wRcorr] = correction(BL, BR, wbar, wL, wR)
if wR < BR
    wRcorr = BR;
    wLcorr = max(BL, 2 * wbar - BR);
elseif wL < BL
    wRcorr = max(BR, 2 * wbar - BL);
    wLcorr = BL;
else
    wRcorr = wR;
    wLcorr = wL;
end