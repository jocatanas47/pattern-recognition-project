clear;
close all;
clc;

rng(100);
N = 500;

M11 = [0 0]';
S11 = [1 0.5; 0.5 1];
M12 = [4 0]';
S12 = [1.2 -0.9; -0.9 1.2];
P11 = 0.75;

M21 = [-2 2]';
S21 = [2 0.4; 0.4 2];
M22 = [2 4]';
S22 = [1.7 -0.5; -0.5 1.7];
P21 = 0.5;

pom = rand(N, 1);

K11 = mvnrnd(M11, S11, N);
K12 = mvnrnd(M12, S12, N);
K1 = (pom < P11).*K11 + (pom >= P11).*K12;

K21 = mvnrnd(M21, S21, N);
K22 = mvnrnd(M22, S22, N);
K2 = (pom < P21).*K21 + (pom >= P21).*K22;

%% dijagram odbiraka
figure();
hold all;
scatter(K1(:, 1), K1(:, 2), 'ro');
scatter(K2(:, 1), K2(:, 2), 'bo');
grid on;
legend("K1", "K2");
title("odbirci");
xlabel("x");
ylabel("y");

%% racunanje fgv
x = -6:0.1:8;
y = -4:0.1:10;
f1 = zeros(length(x), length(y));
f2 = zeros(length(x), length(y));
const11 = 1/(2*pi*det(S11)^0.5);
const12 = 1/(2*pi*det(S12)^0.5);
const21 = 1/(2*pi*det(S21)^0.5);
const22 = 1/(2*pi*det(S22)^0.5);
for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i) y(j)]';
        f11 = const11*exp(-0.5*(X - M11)'*S11^-1*(X - M11));
        f12 = const12*exp(-0.5*(X - M12)'*S12^-1*(X - M12));
        f21 = const21*exp(-0.5*(X - M21)'*S21^-1*(X - M21));
        f22 = const22*exp(-0.5*(X - M22)'*S22^-1*(X - M22));
        f1(i, j) = P11*f11 + (1 - P11)*f12;
        f2(i, j) = P21*f21 + (1 - P21)*f22;
    end
end

%% prikaz histograma i fgv
figure();
subplot(2, 2, 1);
histogram2(K1(:, 1), K1(:, 2), "Normalization", "probability");
xlabel("x");
ylabel("y");
title("histogram K1");
subplot(2, 2, 2);
surf(y, x, f1);
xlabel("x");
ylabel("y");
title("fgv K1");
subplot(2, 2, 3);
histogram2(K2(:, 1), K2(:, 2), "Normalization", "probability");
xlabel("x");
ylabel("y");
title("histogram K2");
subplot(2, 2, 4);
surf(y, x, f2);
xlabel("x");
ylabel("y");
title("fgv K2");

%% Bayes-ov klasifikator minimalne greske
h = -log(f1./f2);

figure();
hold all;
scatter(K1(:, 1), K1(:, 2), 'ro');
scatter(K2(:, 1), K2(:, 2), 'bo');
grid on;
contour(x, y, h', [0 0], 'g', 'Linewidth', 1.5);
xlabel("x");
ylabel("y");
legend("K1", "K2");
title("Bayes-ov klasifikator minimalne greske");

%% greska klasifikacije
greska1 = 0;
for i = 1:length(K1)
    X = K1(i, :)';
    f11 = const11*exp(-0.5*(X - M11)'*S11^-1*(X - M11));
    f12 = const12*exp(-0.5*(X - M12)'*S12^-1*(X - M12));
    f21 = const21*exp(-0.5*(X - M21)'*S21^-1*(X - M21));
    f22 = const22*exp(-0.5*(X - M22)'*S22^-1*(X - M22));
    f1_tren = P11*f11 + (1 - P11)*f12;
    f2_tren = P21*f21 + (1 - P21)*f22;
    if (f2_tren > f1_tren)
        greska1 = greska1 + 1;
    end
end
greska2 = 0;
for i = 1:length(K2)
    X = K2(i, :)';
    f11 = const11*exp(-0.5*(X - M11)'*S11^-1*(X - M11));
    f12 = const12*exp(-0.5*(X - M12)'*S12^-1*(X - M12));
    f21 = const21*exp(-0.5*(X - M21)'*S21^-1*(X - M21));
    f22 = const22*exp(-0.5*(X - M22)'*S22^-1*(X - M22));
    f1_tren = P11*f11 + (1 - P11)*f12;
    f2_tren = P21*f21 + (1 - P21)*f22;
    if (f1_tren > f2_tren)
        greska2 = greska2 + 1;
    end
end

e1_eks = greska1/N
e2_eks = greska2/N

e1_teo = sum(f1(f2 > f1)*0.1*0.1, "all")
e2_teo = sum(f2(f1 > f2)*0.1*0.1, "all")

%% Bayes-ov klasifikator minimalne cene

% c11 = c22 = 0
c12 = 1;
c21 = 10;
h_min_cena = h + log(c12/c21);

figure();
hold all;
scatter(K1(:, 1), K1(:, 2), 'ro');
scatter(K2(:, 1), K2(:, 2), 'bo');
grid on;
contour(x, y, h_min_cena', [0 0], 'g', 'Linewidth', 1.5);
legend("K1", "K2");
xlabel("x");
ylabel("y");
title("Bayes-ov klasifikator minimalne cene sa parametrima c12 = 1 i c21 = 10");

%% Neyman-Pearson-ov klasifikator

% izbor eps0 - greska klasifikacije odbiraka iz druge klase je kod
% Bayes-ovog klasifikatora bila 0.086
% da bi smanjili gresku pri klasifikaciji odbiraka iz prve klase treba da
% povecamo dozvoljenu gresku pri klasifikaciji odbiraka iz druge klase pa
% je
eps0 = 0.12;

mu = 0.01:0.01:10;
eps2 = zeros(1, length(mu));
for br = 1:length(mu)
    mu_tren = mu(br);
    test = h < mu_tren;
    eps2(br) = 0.1*0.1*(sum(f2(h < -log(mu_tren)), "all"));
end

figure();
hold on;
plot(mu, eps2);
plot(mu, ones(1, length(mu))*eps0, '--', 'Linewidth', 1.3);
grid on;
xlabel('mu');
ylabel('eps_{0}');
xlim([0.01, 5]);

ind = find(eps2 < eps0, 1);
e1 = eps2(ind - 1);
e2 = eps2(ind);
mu1 = mu(ind - 1);
mu2 = mu(ind);
% (y - y0) = (y1 - y0)/(x1 - x0)*(x - x0) - interpolacija
mu0 = (eps0 - e1)*(mu2 - mu1)/(e2 - e1) + mu1;

h_np = h + log(mu0);

figure();
hold all;
scatter(K1(:, 1), K1(:, 2), 'ro');
scatter(K2(:, 1), K2(:, 2), 'bo');
grid on;
contour(x, y, h_np', [0 0], 'g', 'Linewidth', 1.5);
legend("K1", "K2");
xlabel("x");
ylabel("y");
title("Neyman-Pearson-ov klasifikator sa usvojenim eps0 = 0.12");

%% Wald-ov sekvencijalni test

% podaci iz prethodnog dela zatatka su bili previse razdvojeni pa se
% Wald-ov test uvek zavrsavao za 1 ili dva koraka -> generisao sam nove
% podatke

rng(100);
N = 500;

M11 = [0 0]';
S11 = [1 0.5; 0.5 1];
M12 = [4 0]';
S12 = [1.2 -0.9; -0.9 1.2];
P11 = 0.75;

M21 = [-2 1]';
S21 = [2 0.4; 0.4 2];
M22 = [2 3]';
S22 = [1.7 -0.5; -0.5 1.7];
P21 = 0.5;

pom = rand(N, 1);

K11 = mvnrnd(M11, S11, N);
K12 = mvnrnd(M12, S12, N);
K1 = (pom < P11).*K11 + (pom >= P11).*K12;

K21 = mvnrnd(M21, S21, N);
K22 = mvnrnd(M22, S22, N);
K2 = (pom < P21).*K21 + (pom >= P21).*K22;

figure();
hold all;
scatter(K1(:, 1), K1(:, 2), 'ro');
scatter(K2(:, 1), K2(:, 2), 'bo');
grid on;
xlabel("x");
ylabel("y");
legend("K1", "K2");
title("odbirci korisceni za Wald-ov test");

e1 = 0.01:0.005:0.08;
e2 = 0.01:0.005:0.08;

potrebni_odbirci = zeros(length(e1), length(e2));
% razmatramo klasu K1

const11 = 1/(2*pi*det(S11)^0.5);
const12 = 1/(2*pi*det(S12)^0.5);
const21 = 1/(2*pi*det(S21)^0.5);
const22 = 1/(2*pi*det(S22)^0.5);

for i = 1:length(e1)
    for j = 1:length(e2)
        A = -(e1(i) - 1)/e2(j);
        B = -e1(i)/(e2(j) - 1);
        lm = 1;
        br = 0;
        for k = 1:size(K1, 1)
            br = br + 1;
            K1_tren = K1(k, :);
            f11 = const11*exp(-0.5*(K1_tren' - M11)'*S11^-1*(K1_tren' - M11));
            f12 = const12*exp(-0.5*(K1_tren' - M12)'*S12^-1*(K1_tren' - M12));
            f21 = const21*exp(-0.5*(K1_tren' - M21)'*S21^-1*(K1_tren' - M21));
            f22 = const22*exp(-0.5*(K1_tren' - M22)'*S22^-1*(K1_tren' - M22));
            f1 = P11*f11 + (1 - P11)*f12;
            f2 = P21*f21 + (1 - P21)*f22;
            lm = lm * f1/f2;
            if (lm >= A)
                potrebni_odbirci(i, j) = br;
                break;
            end
        end
    end
end

figure();
surf(e2, e1, potrebni_odbirci);
xlabel("e2");
ylabel("e1");
title("potreban broj odbiraka");