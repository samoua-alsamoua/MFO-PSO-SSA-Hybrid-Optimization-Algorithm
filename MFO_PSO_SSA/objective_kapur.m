function [ fitness ] = objective_kapur(p, x, level)
    rgbValue = 256;
    fitness = 0;
    pArray = p(1 : x(1));
    fitness = fitness + calculate(pArray);
    for j = 2 : level
        pArray = p(x(j-1) + 1 : x(j));
        fitness = fitness + calculate(pArray);
    end
    pArray = p(x(level) + 1 : rgbValue);
    fitness = fitness + calculate(pArray);
    if fitness == 0
        fitness = 10E-100;
    end
end

function fValue = calculate(pArray)
    index = pArray == 0;
    index = index .* eps;
    pArray = pArray + index;
    sumP = sum(pArray);
    pSum = pArray / sumP;
    fValue = -sum(pSum .* log2(pSum));
end