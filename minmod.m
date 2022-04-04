function S = minmod(A,B,C)
S = min(min(A,B),C) .* and(and(A>0,B>0),C>0) + max(max(A,B),C) .* and(and(A<0,B<0),C<0);
end