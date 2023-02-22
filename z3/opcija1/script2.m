clear;
close all;
clc;

N = 500;

phi1 = 2*pi*rand(1, N);
phi1 = phi1.^2;
rho1 = rand(1, N);
X = zeros(2, N);
X(1, :) = rho1.*cos(phi1);
X(2, :) = rho1.*sin(phi1);

pom = rand(1, N);
phi2 = pi/4*rand(1, N).*(pom < 0.5) + (pi + pi/4*rand(1, N)).*(pom >= 0.5);
rho2 = rand(1, N) + 1.3;
Y = zeros(2, N);
Y(1, :) = rho2.*cos(phi2);
Y(2, :) = rho2.*sin(phi2);

figure();
hold on;
scatter(X(1, :), X(2, :), 'bo');
scatter(Y(1, :), Y(2, :), 'ro');
xlabel("x");
ylabel("y");

Gamma = [ones(N, 1); 2*ones(N, 1)];

U = [-1*ones(1, N), ones(1, N); ...
    -1*X, Y; ...
    -1*X(1, :).^2, Y(1, :).^2; ...
    -1*X(2, :).^2, Y(2, :).^2; ...
    -2*X(1, :).*X(2, :), 2*Y(1, :).*Y(2, :)];
W = (U*U')^(-1)*U*Gamma;
v0 = W(1);
V1 = W(2);
V2 = W(3);
Q11 = W(4);
Q22 = W(5);
Q12 = W(6);

x1 = -3:0.1:3;
x2 = -3:0.1:3;
[x1, x2] = meshgrid(x1, x2);
h = v0 + V1*x1 + V2*x2 + Q11*x1.^2 + Q22*x2.^2 + Q12*x1.*x2;

contour(x1, x2, h, [0, 0]);
legend("K1", "K2", "klasifikaciona linija");
title("kvadratni klasifikator metodom zeljenog izlaza");