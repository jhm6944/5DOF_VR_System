function blended_img = blend(im0, im1)
height = size(im0, 1);
width = size(im0, 2);

map0 = repmat((width-1:-1:0)/(width-1), [height, 1]);
map1 = repmat((0:width-1)/(width-1), [height, 1]);

blended_img = im2uint8(map0 .* im2double(im0) + map1 .* im2double(im1));
end