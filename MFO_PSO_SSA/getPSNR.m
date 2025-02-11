function [psnr] = getPSNR(image, test)
    distance = double(image - test) .^ 2;
    sumDist = sum(sum(distance)); 
    if( sumDist <= 1e-10 ) 
        psnr = 0;
    else
        mse  = sumDist / double(size(image, 1) * size(image, 2) * size(image, 3));
        psnr = 10.0 * log10((255 * 255) / mse);
    end
end