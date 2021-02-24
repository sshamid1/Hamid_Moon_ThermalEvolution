function [fchi] = fchi(x,ri,Lp)
%FCHI "Yet another useful function"

fchi = x^3 * (-1/3*(ri/Lp)^2 + 0.2*(1+(ri/Lp)^2)*x^2 - 13/70*x^4);

end

