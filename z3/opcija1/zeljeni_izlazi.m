function [V, v0] = zeljeni_izlazi(X1, X2, N, Gamma)
U = [-ones(1, N), ones(1, N); ...
    -1*X1, X2];
W = (U*U')^(-1)*U*Gamma;
v0 = W(1);
V = [W(2); W(3)];
end