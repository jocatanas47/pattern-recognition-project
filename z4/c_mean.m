function [klasteri, broj_iteracija] = c_mean(X, br_klastera)
klasteri = randi([1, br_klastera], 1, size(X, 2));
centri = zeros(2, br_klastera);
for i = 1:br_klastera
    centri(:, i) = mean(X(:, klasteri == i), 2);
end
stari_klasteri = klasteri;

broj_iteracija = 0;
while (true)
    broj_iteracija = broj_iteracija + 1;
    razdaljine = zeros(br_klastera, size(X, 2));
    for i = 1:br_klastera
        for j = 1:size(X, 2)
            razdaljine(i, j) = sum((X(:, j) - centri(:, i)).^2, 1);
        end
    end
    [~, klasteri] = min(razdaljine, [], 1);
    if (isequal(stari_klasteri, klasteri))
        break;
    end
    stari_klasteri = klasteri;
    for i = 1:br_klastera
        centri(:, i) = mean(X(:, klasteri == i), 2);
    end
end
end