clear;
close all;
clc;

%% generisanje klasa

N = 500;

M1 = [0; 0];
S1 = [1, 0; 0, 1];
M2 = [6; 6];
S2 = [1.2, 0.7; 0.7, 1.2];
M3 = [-3; 4];
S3 = [1, 0.2; 0.2, 1];
M4 = [-1; -6];
S4 = [1, 0.3; 0.3, 1];

X1 = mvnrnd(M1, S1, N)';
X2 = mvnrnd(M2, S2, N)';
X3 = mvnrnd(M3, S3, N)';
X4 = mvnrnd(M4, S4, N)';

figure();
hold all;
scatter(X1(1, :), X1(2, :));
scatter(X2(1, :), X2(2, :));
scatter(X3(1, :), X3(2, :));
scatter(X4(1, :), X4(2, :));
legend("K1", "K2", "K3", "K4");
xlabel("x");
ylabel("y");
title("pocetne klase");

%% c-mean klasterizacija
X = [X1, X2, X3, X4];
br_iteracija = 10;
trajanja = zeros(1, br_iteracija);

for k = 1:br_iteracija
    [klasteri, trajanja(k)] = c_mean(X, 4);
end

figure();
plot(trajanja);
ylabel("potreban broj iteracija");
title("oetljivost na pocetnu klasterizaciju c-mean metode");
sr_trajanje_c_mean = mean(trajanja)

K1 = X(:, klasteri == 1);
K2 = X(:, klasteri == 2);
K3 = X(:, klasteri == 3);
K4 = X(:, klasteri == 4);
figure();
hold all;
scatter(K1(1, :), K1(2, :));
scatter(K2(1, :), K2(2, :));
scatter(K3(1, :), K3(2, :));
scatter(K4(1, :), K4(2, :));
legend("K1", "K2", "K3", "K4");
xlabel("x");
ylabel("y");
title("rezultat c-mean klasterizacije");

%% nepoznat broj klastera
max_broj_klastera = 6;
J = zeros(1, max_broj_klastera);

for i = 1:max_broj_klastera
    [klasteri, ~] = c_mean(X, i);
    P_est = zeros(1, i);
    M_est = zeros(2, i);
    S_est = zeros(2, 2, i);
    for j = 1:i
        P_est(j) = size(X(:, klasteri == j), 2)/size(X, 2);
        M_est(:, j) = mean(X(:, klasteri(:) == j), 2);
        S_est(:, :, j) = cov(X(:, klasteri == j)');
    end
    M0 = zeros(2, 1);
    for j = 1:i
        M0 = M0 + P_est(j)*M_est(:, j);
    end
    Sb = zeros(2, 2);
    Sw = zeros(2, 2);
    for j = 1:i
        Sb = Sb + P_est(j)*(M_est(:, j) - M0)*(M_est(:, j) - M0)';
        Sw = Sw + P_est(j)*S_est(:, :, j);
    end
    Sm = Sb + Sw;
    J(i) = trace(Sm^-1*Sw);
end
    
figure();
plot(1:max_broj_klastera, J);
xlabel("broj klastera");
ylabel("J");
title("vrednost J za razlicite brojeve klastera, c-mean");
grid on;

%% metod maksimalne verodostojnosti
X = [X1, X2, X3, X4];
br_iteracija = 10;
trajanja = zeros(1, br_iteracija);

for k = 1:br_iteracija
    [klasteri, trajanja(k)] = maksimalna_verodostojnost(X, 4);
end

figure();
plot(trajanja);
ylabel("potreban broj iteracija");
sr_trajanje_maks_ver = mean(trajanja)
title("osetljivost ml metodom");

K1 = X(:, klasteri == 1);
K2 = X(:, klasteri == 2);
K3 = X(:, klasteri == 3);
K4 = X(:, klasteri == 4);
figure();
hold all;
scatter(K1(1, :), K1(2, :));
scatter(K2(1, :), K2(2, :));
scatter(K3(1, :), K3(2, :));
scatter(K4(1, :), K4(2, :));
legend("K1", "K2", "K3", "K4");
xlabel("x");
ylabel("y");
title("rezultat klasterizacije ml metodom");

%% nepoznat broj klastera
max_broj_klastera = 6;
J = zeros(1, max_broj_klastera);

for i = 1:max_broj_klastera
    [klasteri, ~] = maksimalna_verodostojnost(X, i);
    P_est = zeros(1, i);
    M_est = zeros(2, i);
    S_est = zeros(2, 2, i);
    for j = 1:i
        P_est(j) = size(X(:, klasteri == j), 2)/size(X, 2);
        M_est(:, j) = mean(X(:, klasteri(:) == j), 2);
        S_est(:, :, j) = cov(X(:, klasteri == j)');
    end
    M0 = zeros(2, 1);
    for j = 1:i
        M0 = M0 + P_est(j)*M_est(:, j);
    end
    Sb = zeros(2, 2);
    Sw = zeros(2, 2);
    for j = 1:i
        Sb = Sb + P_est(j)*(M_est(:, j) - M0)*(M_est(:, j) - M0)';
        Sw = Sw + P_est(j)*S_est(:, :, j);
    end
    Sm = Sb + Sw;
    J(i) = trace(Sm^-1*Sw);
end
    
figure();
plot(1:max_broj_klastera, J);
xlabel("broj klastera");
ylabel("J");
title("vrednost J za razlicite brojeve klastera, ml metoda");
grid on;

%% klasterizacija nelinearno separabilnih odbiraka metodom maksimalne verodostojnosti
N = 500;

phi1 = 2*pi*rand(1, N);
phi1 = phi1.^2;
rho1 = rand(1, N);
X = zeros(2, N);
X(1, :) = rho1.*cos(phi1);
X(2, :) = rho1.*sin(phi1);

pom = rand(1, N);
phi2 = 2*pi*rand(1, N);
rho2 = rand(1, N) + 2.5;
Y = zeros(2, N);
Y(1, :) = rho2.*cos(phi2);
Y(2, :) = rho2.*sin(phi2);

figure();
hold all;
scatter(X(1, :), X(2, :));
scatter(Y(1, :), Y(2, :));
legend("K1", "K2");
xlabel("x");
ylabel("y");
title("pocetne nelinearno separabilne klase");

podaci = [X, Y];

br_iteracija = 10;
trajanja = zeros(1, br_iteracija);
for k = 1:br_iteracija
    [klasteri, trajanja(k)] = maksimalna_verodostojnost(podaci, 2);
end

figure();
plot(trajanja);
ylabel("potreban broj iteracija");
sr_trajanje_maks_ver_nelinearno_separabilne = mean(trajanja)
title("osetljivost, ml metoda, nelinearno sep klase");

K1 = podaci(:, klasteri == 1);
K2 = podaci(:, klasteri == 2);
figure();
hold all;
scatter(K1(1, :), K1(2, :));
scatter(K2(1, :), K2(2, :));
legend("K1", "K2");
xlabel("x");
ylabel("y");
title("klasterizacija nelinearno separabilnih odbiraka ml metodom");

%% nepoznat broj klastera
max_broj_klastera = 3;
J = zeros(1, max_broj_klastera);

for i = 1:max_broj_klastera
    [klasteri, ~] = maksimalna_verodostojnost(podaci, i);
    P_est = zeros(1, i);
    M_est = zeros(2, i);
    S_est = zeros(2, 2, i);
    for j = 1:i
        P_est(j) = size(podaci(:, klasteri == j), 2)/size(podaci, 2);
        M_est(:, j) = mean(podaci(:, klasteri(:) == j), 2);
        S_est(:, :, j) = cov(podaci(:, klasteri == j)');
    end
    M0 = zeros(2, 1);
    for j = 1:i
        M0 = M0 + P_est(j)*M_est(:, j);
    end
    Sb = zeros(2, 2);
    Sw = zeros(2, 2);
    for j = 1:i
        Sb = Sb + P_est(j)*(M_est(:, j) - M0)*(M_est(:, j) - M0)';
        Sw = Sw + P_est(j)*S_est(:, :, j);
    end
    Sm = Sb + Sw;
    J(i) = trace(Sm^-1*Sw);
end
    
figure();
plot(1:max_broj_klastera, J);
xlabel("broj klastera");
ylabel("J");
title("vrednost J za razlicite brojeve klastera, ml metoda, nelinearno sep klase");
grid on;