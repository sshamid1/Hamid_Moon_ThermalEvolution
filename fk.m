function [fk] = fk(x,Ap)
%FK "A Useful Function"

fk = 0.2*x^5 * (1 + 5/7*(2+4*Ap)*x^2 + 5/9 * (3 + 10*Ap + 4*Ap^2)*x^4);

end

