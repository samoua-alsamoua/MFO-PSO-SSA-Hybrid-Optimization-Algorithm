function [ fitness ] = testFunction( x, p )

global count;
[~, counter] = size(x);
count = count + counter;

% OTSU METHOD
%     [dimension, pop] = size(x);
%     fitness = zeros(1, pop);
%     for i = 1 : pop
%         fitness(i) = (1 / objective_otsu(p, fix(sort(x(:,i))), dimension));
%     end

% KAPUR METHOD    
    [dimension, pop] = size(x);
    fitness = zeros(1, pop);
    for i = 1 : pop
        fitness(i) = (1 / objective_kapur(p, fix(sort(x(:,i))), dimension));
    end

end

