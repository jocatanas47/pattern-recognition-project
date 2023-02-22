function obelezje = obelezje2(img)
a = sum(img(:, 1:round(end/2)), "all");
b = sum(img(:, round(end/2):end), "all");
obelezje = a/b;
end