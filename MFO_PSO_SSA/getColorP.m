function [ p ] = getColorP(colorImage, index)

    image = colorImage(:, :, index);
    [row, column] = size(image);
    pixel = row * column;
    vector = im2uint8(image(:));
    histogram = imhist(vector, 256);
    p = (histogram / pixel)';

end