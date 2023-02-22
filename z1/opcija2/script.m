clear;
close all;
clc;

%% makaze
list = dir("scissors_predobradjen/*.png");
filenames = string({list.name});

n = round(0.8*length(filenames));

obelezja_train_s = zeros(n, 2);
for i = 1:n
    img = imread("scissors_predobradjen/" + filenames(i));
    obelezja_train_s(i, 1) = obelezje1(img);
    obelezja_train_s(i, 2) = obelezje2(img);
    %obelezja_train_s(i, 3) = obelezje3(img);
end

obelezja_test_s = zeros(length(filenames) - n, 2);
for i = 1:length(filenames) - n
    img = imread("scissors_predobradjen/" + filenames(i + n));
    obelezja_test_s(i, 1) = obelezje1(img);
    obelezja_test_s(i, 2) = obelezje2(img);
    %obelezja_test_s(i, 3) = obelezje3(img);
end

%% kamen
list = dir("rock_predobradjen/*.png");
filenames = string({list.name});

n = round(0.8*length(filenames));

obelezja_train_r = zeros(n, 2);
for i = 1:n
    img = imread("rock_predobradjen/" + filenames(i));
    obelezja_train_r(i, 1) = obelezje1(img);
    obelezja_train_r(i, 2) = obelezje2(img);
    %obelezja_train_r(i, 3) = obelezje3(img);
end

obelezja_test_r = zeros(length(filenames) - n, 2);
for i = 1:length(filenames) - n
    img = imread("rock_predobradjen/" + filenames(i + n));
    obelezja_test_r(i, 1) = obelezje1(img);
    obelezja_test_r(i, 2) = obelezje2(img);
    %obelezja_test_r(i, 3) = obelezje3(img);
end

%% papir

list = dir("paper_predobradjen/*.png");
filenames = string({list.name});

n = round(0.8*length(filenames));

obelezja_train_p = zeros(n, 2);
for i = 1:n
    img = imread("paper_predobradjen/" + filenames(i));
    obelezja_train_p(i, 1) = obelezje1(img);
    obelezja_train_p(i, 2) = obelezje2(img);
    %obelezja_train_p(i, 3) = obelezje3(img);
end

obelezja_test_p = zeros(length(filenames) - n, 2);
for i = 1:length(filenames) - n
    img = imread("paper_predobradjen/" + filenames(i + n));
    obelezja_test_p(i, 1) = obelezje1(img);
    obelezja_test_p(i, 2) = obelezje2(img);
    %obelezja_test_p(i, 3) = obelezje3(img);
end

%% klasifikacija zasnovana na testiranju hipoteza

P_s = size(obelezja_train_s, 1);
P_r = size(obelezja_train_r, 1);
P_p = size(obelezja_train_p, 1);
k = P_s + P_r + P_p;
P_s = P_s/k;
P_r = P_r/k;
P_p = P_p/k;
[f_s, ivicex_s, ivicey_s] = histcounts2(obelezja_train_s(:, 1), obelezja_train_s(:, 2), 25, "Normalization", "probability");
[f_r, ivicex_r, ivicey_r] = histcounts2(obelezja_train_r(:, 1), obelezja_train_r(:, 2), 25, "Normalization", "probability");
[f_p, ivicex_p, ivicey_p] = histcounts2(obelezja_train_p(:, 1), obelezja_train_p(:, 2), 25, "Normalization", "probability");

X_test = [obelezja_test_s; obelezja_test_r; obelezja_test_p];
Y_test = [ones(1, size(obelezja_test_s, 1)), 2*ones(1, size(obelezja_test_r, 1)), 3*ones(1, size(obelezja_test_p, 1))];
pred = zeros(1, length(Y_test));
for i = 1:length(pred)
    xbin = find(ivicex_s <= X_test(i, 1), 1, "last");
    ybin = find(ivicey_s <= X_test(i, 2), 1, "last");
    if (xbin > 25)
        xbin = 25;
    end
    if (ybin > 25)
        ybin = 25;
    end
    val_s = f_s(xbin, ybin);
    
    xbin = find(ivicex_r <= X_test(i, 1), 1, "last");
    ybin = find(ivicey_r <= X_test(i, 2), 1, "last");
    if (xbin > 25)
        xbin = 25;
    end
    if (ybin > 25)
        ybin = 25;
    end
    val_r = f_r(xbin, ybin);
    
    xbin = find(ivicex_p <= X_test(i, 1), 1, "last");
    ybin = find(ivicey_p <= X_test(i, 2), 1, "last");
    if (xbin > 25)
        xbin = 25;
    end
    if (ybin > 25)
        ybin = 25;
    end
    val_p = f_p(xbin, ybin);
    
    vector = [P_s*val_s, P_r*val_r, P_p*val_p];
    [~, pred(i)] = max(vector);
end

cm1 = confusionmat(Y_test, pred)

%% histogrami slova
figure()
hist3(obelezja_train_s);
title("histogram obelezja Makaza");

figure()
hist3(obelezja_train_r);
title("histogram obelezja Kamena");

%% crtanje i parametarski klasifikator
figure();
hold on;
scatter(obelezja_train_s(:, 1), obelezja_train_s(:, 2));
scatter(obelezja_train_r(:, 1), obelezja_train_r(:, 2));
%scatter(obelezja_train_p(:, 1), obelezja_train_p(:, 2));
xlabel("x");
ylabel("y");
ylim([0, 15]);

Ms_est = mean(obelezja_train_s, 1)';
Mr_est = mean(obelezja_train_r, 1)';
Ss_est = cov(obelezja_train_s);
Sr_est = cov(obelezja_train_r);

[V, v0_opt] = linearni_klasifikator(obelezja_train_s', obelezja_train_r', ...
    Ms_est, Mr_est, Ss_est, Sr_est);

x = 0:0.01:1.2;
y = -(v0_opt + V(1)*x)/V(2);

plot(x, y, "Linewidth", 1.4);

legend("scissors", "rock", "klasifikaciona prava");

X_test = [obelezja_test_s; obelezja_test_r];
Y_test = [-ones(1, size(obelezja_test_s, 1)), ones(1, size(obelezja_test_r, 1))];
pred = zeros(1, length(Y_test));

for i = 1:length(pred)
    pred(i) = klasifikuj(X_test(i, :)', V, v0_opt);
end

cm2 = confusionmat(Y_test, pred)