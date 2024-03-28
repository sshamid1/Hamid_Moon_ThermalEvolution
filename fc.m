function [fc] = fc(x, delta, Ap)
%FC "The last useful function, for now"

fc = x^3 * (1 - 0.6*(delta+1)*x^2 - 3/14*(delta+1)*(2*Ap-delta)*x^4);

end

