function obelezje = obelezje3(img)
% ruka = bwconvhull(img);
% stats = regionprops(ruka, 'Centroid');
% centar = stats.Centroid;
% [rows, cols] = find(ruka);
% distances = sqrt((rows - centar(2)).^2 + (cols - centar(1)).^2);
% obelezje = mean(distances);
stats = regionprops(img, "Area");
obelezje = stats.Area;
end