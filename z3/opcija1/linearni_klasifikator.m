function [V, v0] = linearni_klasifikator(X1, X2, M1, M2, S1, S2, N)
s = 0:1e-3:1;
v0_opt_s = zeros(1, length(s));
Neps_s = zeros(1, length(s));
for i = 1:length(s)
    V = ((s(i)*S1 + (1 - s(i))*S2)^(-1))*(M2 - M1);
    Y1 = V'*X1;
    Y2 = V'*X2;
    Y = [Y1, Y2];
    Y = sort(Y);
    v0 = zeros(1, length(Y) - 1);
    Neps = zeros(1, length(Y) - 1);
    for j = 1:(length(Y) - 1)
        v0(j) = -(Y(j) + Y(j + 1))/2;
        Neps(j) = 0;
        for k = 1:N
            if (Y1(k) > -v0(j))
                Neps(j) = Neps(j) + 1;
            end
        end
        for k = 1:N
            if (Y2(k) < -v0(j))
                Neps(j) = Neps(j) + 1;
            end
        end
    end
    [Neps_s(i), indx] = min(Neps);
    v0_opt_s(i) = v0(indx);
end

[~, indx] = min(Neps_s);
v0 = v0_opt_s(indx);
s_opt = s(indx);
V = ((s_opt*S1 + (1 - s_opt)*S2)^-1)*(M2 - M1);
end