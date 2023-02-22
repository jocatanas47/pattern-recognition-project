clear;
close all;
clc;

%% generisanje klasa
N = 500;
rng(100);

M1 = [1; -1];
S1 = [1, 0; 0, 1];
M2 = [5; 8];
S2 = [1.2, 0.4; 0.4, 1.2];
M3 = [9; 3];
S3 = [1, 0.3; 0.3, 1];

X1 = mvnrnd(M1, S1, N)';
X2 = mvnrnd(M2, S2, N)';
X3 = mvnrnd(M3, S3, N)';

figure();
hold all;
scatter(X1(1, :), X1(2, :), 'ro');
scatter(X2(1, :), X2(2, :), 'bo');
scatter(X3(1, :), X3(2, :), 'go');
legend("K1", "K2", "K3");
xlabel("x");
ylabel("y");
title("generisani odbirci");

%% linearni klasifikator - druga numericka metoda
M1_est = mean(X1, 2);
M2_est = mean(X2, 2);
M3_est = mean(X3, 2);
S1_est = cov(X1');
S2_est = cov(X2');
S3_est = cov(X3');

[V12, v0_opt12] = linearni_klasifikator(X1, X2, M1_est, M2_est, S1_est, S2_est, N); 
[V23, v0_opt23] = linearni_klasifikator(X2, X3, M2_est, M3_est, S2_est, S3_est, N);
[V31, v0_opt31] = linearni_klasifikator(X3, X1, M3_est, M1_est, S3_est, S1_est, N);

x = -2:0.1:12;
x12 = -(v0_opt12 + V12(1)*x)/V12(2);
x23 = -(v0_opt23 + V23(1)*x)/V23(2);
x31 = -(v0_opt31 + V31(1)*x)/V31(2);

figure();
hold all;
scatter(X1(1, :), X1(2, :), 'ro');
scatter(X2(1, :), X2(2, :), 'bo');
scatter(X3(1, :), X3(2, :), 'go');
plot(x, x12);
plot(x, x23);
plot(x, x31);
xlabel("x");
ylabel("y");
legend("K1", "K2", "K3", "x12", "x23", "x31");
title("linearna klasifikacija - druga numericka metoda");

Y_true = [ones(1, N), 2*ones(1, N), 3*ones(1, N)];
Y_pred = zeros(1, length(Y_true));
k = zeros(3, length(Y_true));
X = [X1, X2, X3];
for i = 1:length(X)
    k12 = klasifikuj(X(:, i), V12, v0_opt12);
    k23 = klasifikuj(X(:, i), V23, v0_opt23);
    k31 = klasifikuj(X(:, i), V31, v0_opt31);
    k(:, i) = [k12; k23; k31];
    if (k31 == 1 && k12 == -1)
        Y_pred(i) = 1;
    end
    if (k12 == 1 && k23 == -1)
        Y_pred(i) = 2;
    end
    if (k23 == 1 && k31 == -1)
        Y_pred(i) = 3;
    end
end

cm_lk = confusionmat(Y_true, Y_pred)

%% metod zeljenih izlaza

Gamma12 = [ones(N, 1); 2*ones(N, 1)];
Gamma23 = [3*ones(N, 1); ones(N, 1)];
Gamma31 = [7*ones(N, 1); ones(N, 1)];

[V12, v0_opt12] = zeljeni_izlazi(X1, X2, N, Gamma12); 
[V23, v0_opt23] = zeljeni_izlazi(X2, X3, N, Gamma23);
[V31, v0_opt31] = zeljeni_izlazi(X3, X1, N, Gamma31);

x = -2:0.1:12;
x12 = -(v0_opt12 + V12(1)*x)/V12(2);
x23 = -(v0_opt23 + V23(1)*x)/V23(2);
x31 = -(v0_opt31 + V31(1)*x)/V31(2);

figure();
hold all;
scatter(X1(1, :), X1(2, :), 'ro');
scatter(X2(1, :), X2(2, :), 'bo');
scatter(X3(1, :), X3(2, :), 'go');
plot(x, x12);
plot(x, x23);
plot(x, x31);
xlabel("x");
ylabel("y");
legend("K1", "K2", "K3", "x12", "x23", "x31");
title("metod zeljenih izlaza");

Y_true = [ones(1, N), 2*ones(1, N), 3*ones(1, N)];
Y_pred = zeros(1, length(Y_true));
k = zeros(3, length(Y_true));
X = [X1, X2, X3];
for i = 1:length(X)
    k12 = klasifikuj(X(:, i), V12, v0_opt12);
    k23 = klasifikuj(X(:, i), V23, v0_opt23);
    k31 = klasifikuj(X(:, i), V31, v0_opt31);
    k(:, i) = [k12; k23; k31];
    if (k31 == 1 && k12 == -1)
        Y_pred(i) = 1;
    end
    if (k12 == 1 && k23 == -1)
        Y_pred(i) = 2;
    end
    if (k23 == 1 && k31 == -1)
        Y_pred(i) = 3;
    end
end

cm_zi = confusionmat(Y_true, Y_pred)