clear;
close all;
clc;

str = "paper";
list = dir(str + "/*.png");
filenames = string({list.name});

for i = 1:length(filenames)
    img = imread(str + "/" + filenames(i));
    
    % binarizacija
    r = img(:, :, 1);
    g = img(:, :, 2);
    b = img(:, :, 3);
    mask = (g > r) & (g > b);
    
    % izdvajanje bitnog dela
    edges = edge(mask, "canny");
    L = logical(edges);
    stats = regionprops(L, "all");
    [~, max_index] = max([stats.Area]);
    useful_area = stats(max_index);
    bounding_box = useful_area.BoundingBox;
    img_crop = imcrop(mask, bounding_box);
    
    % cuvanje slike
    imwrite(img_crop, str + "_predobradjen/" + filenames(i));
end