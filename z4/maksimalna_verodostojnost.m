function [klasteri, broj_iteracija] = maksimalna_verodostojnost(X, br_klastera)
T = 10^-3;

klasteri = randi([1, br_klastera], 1, size(X, 2));

P_l = zeros(1, br_klastera);
M_l = zeros(2, br_klastera);
S_l = zeros(2, 2, br_klastera);
q_l = zeros(size(X, 2), br_klastera);
q_l1 = zeros(size(X, 2), br_klastera);
f_k = zeros(size(X, 2), br_klastera);
f = zeros(size(X, 2), 1);

for i = 1:br_klastera
    P_l(i) = sum(klasteri(:) == i, 1)/size(X, 2);
    M_l(:, i) = mean(X(:, klasteri == 1), 2);
    S_l(:, :, i) = cov(X(:, klasteri(:) == i)');
end

for i = 1:br_klastera
    k1 = 1/(2*pi*det(S_l(:, :, i))^0.5);
    k2 = S_l(:, :, i)^-1;
    for j = 1:size(X, 2)
        f_k(j, i) = k1 * exp(-0.5*(X(:, j) - M_l(:, i))'*k2*(X(:, j) - M_l(:, i)));
    end
end
for i = 1:br_klastera
    f = f + P_l(i)*f_k(:, i);
end
for i = 1:br_klastera
    q_l(:, i) = P_l(i) * f_k(:, i)./f;
end

broj_iteracija = 0;
while (broj_iteracija < 300)
    broj_iteracija = broj_iteracija + 1;
    
    S_l = zeros(2, 2, br_klastera);
    P_l = sum(q_l, 1)/size(X, 2);
    for i = 1:br_klastera
        M_l(:, i) = (q_l(:, i)'*X')'/size(X, 2)/P_l(i);
        for j = 1:size(X, 2)
            S_l(:, :, i) = S_l(:, :, i) + ...
                q_l(j, i)*(X(:, j) - M_l(:, i))*(X(:, j) - M_l(:, i))';
        end
        S_l(:, :, i) = S_l(:, :, i)/size(X, 2)/P_l(i);
    end
    
    for i = 1:br_klastera
        k1 = 1/(2*pi*det(S_l(:, :, i))^0.5);
        k2 = S_l(:, :, i)^-1;
        for j = 1:size(X, 2)
            f_k(j, i) = k1 * exp(-0.5*(X(:, j) - M_l(:, i))'*k2*(X(:, j) - M_l(:, i)));
        end
    end
    
    f = zeros(size(X, 2), 1);
    for i = 1:br_klastera
        f = f + P_l(i)*f_k(:, i);
    end
    for i = 1:br_klastera
        q_l1(:, i) = P_l(i) .* f_k(:, i)./f;
    end
    
    razlika = abs(q_l - q_l1);
    if (max(razlika(:)) < T)
        break;
    end
    
    q_l = q_l1;
end
[~, klasteri] = max(q_l');
end