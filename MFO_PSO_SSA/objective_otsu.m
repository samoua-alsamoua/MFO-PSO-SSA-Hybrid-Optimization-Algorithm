function [ variance ] = objective_otsu(p, x, level)
    rgbValue = 256;
    arrayX = (1:x(1));
    arrayMax = (1:rgbValue);
    maxMultiple = arrayMax .* p(arrayMax);
    sumMax = sum(maxMultiple);
    pArray = p(arrayX);
    sumP = sum(pArray);
    if sumP == 0
        variance = 0;
    else
        variance = sumP * (sum(arrayX .* pArray / sumP) - sumMax) ^ 2; 
    end
    for i = 2 : level 
        arrayXi = x(i-1) + 1 : x(i);
        piArray = p(arrayXi);
        sumPi = sum(piArray);
        if sumPi ~= 0
            variance = variance + sumPi * (sum(arrayXi .* piArray / sumPi) - sumMax) ^ 2;
        end        
    end
    arrayXl = x(level) + 1 : rgbValue;
    plArray = p(arrayXl);
    sumPl = sum(plArray);
    if sumPl ~= 0
        variance = variance + sumPl * (sum(arrayXl .* plArray / sumPl) - sumMax) ^ 2;
    end
    if variance == 0
        variance = 10E-100;
    end
end