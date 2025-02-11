function [bestSolution, bestFitness, iter] = AVOA(p, dimension, maxIteration)

    settingAlg;
    Max_iter = maxIteration;
    dim = dimension;
    lb = lbArray;
    ub = ubArray;
    N = 30;        %%%% pop_size
    
    % initialize Best_vulture1, Best_vulture2
    Best_vulture1_X = zeros(1, dim);
    Best_vulture1_F = inf;
    Best_vulture2_X = zeros(1, dim);
    Best_vulture2_F = inf;

    % Initialize the first random population of vultures
    Positions = initialization(N, dim, ub, lb);

    %% Controlling parameter
    p1 = 0.6;
    p2 = 0.4;
    p3 = 0.6;
    alpha = 0.8;
    betha = 0.2;
    gamma = 2.5;

    %% Main loop
    l = 1; % Loop counter

    while l < Max_iter
        for i = 1:size(Positions, 1)
            % Calculate the fitness of the population
            PositionsFitness(1, i) = testFunction(Positions(i, :)', p);
            l = l + 1;
            
            % Update the first best two vultures if needed
            if PositionsFitness(1, i) < Best_vulture1_F
                Best_vulture1_F = PositionsFitness(1, i); % Update the first best vulture
                Best_vulture1_X = Positions(i, :);
            end
            if PositionsFitness(1, i) > Best_vulture1_F && PositionsFitness(1, i) < Best_vulture2_F
                Best_vulture2_F = PositionsFitness(1, i); % Update the second best vulture
                Best_vulture2_X = Positions(i, :);
            end
        end

        a = unifrnd(-2, 2, 1, 1) * ((sin((pi / 2) * (l / Max_iter)) ^ gamma) + cos((pi / 2) * (l / Max_iter)) - 1);
        P1 = (2 * rand + 1) * (1 - (l / Max_iter)) + a;

        % Update the location
        for i = 1:size(Positions, 1)
            Positions(i, :) = Positions(i, :); % pick the current vulture back to the population
            F = P1 * (2 * rand() - 1);

            random_vulture_X = random_select(Best_vulture1_X, Best_vulture2_X, alpha, betha);

            if abs(F) >= 1 % Exploration:
                Positions(i, :) = exploration(Positions(i, :), random_vulture_X, F, p1, ub, lb);
            elseif abs(F) < 1 % Exploitation:
                Positions(i, :) = exploitation(Positions(i, :), Best_vulture1_X, Best_vulture2_X, random_vulture_X, F, p2, p3, dim, ub, lb);
            end

            Positions(i, :) = Positions(i, :); % place the current vulture back into the population
        end

        Positions = boundaryCheck(Positions, lb, ub);

    end

    bestSolution = Best_vulture1_X;
    bestFitness = Best_vulture1_F;
    iter = l;
end

% 
% % This function initialize the first population of search agents
function [ X ]=initialization(N,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(N,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end

function [random_vulture_X]=random_select(Best_vulture1_X,Best_vulture2_X,alpha,betha)

    probabilities=[alpha, betha ];
    
    if (rouletteWheelSelection( probabilities ) == 1)
            random_vulture_X=Best_vulture1_X;
    else
            random_vulture_X=Best_vulture2_X;
    end

end

function [index] = rouletteWheelSelection(x)

    index=find(rand() <= cumsum(x) ,1,'first');

end

function [current_vulture_X] = exploitation(current_vulture_X, Best_vulture1_X, Best_vulture2_X, ...
                                                                      random_vulture_X, F, p2, p3, variables_no, upper_bound, lower_bound)

% phase 1
    if  abs(F)<0.5
        if rand<p2
            A=Best_vulture1_X-((Best_vulture1_X.*current_vulture_X)./(Best_vulture1_X-current_vulture_X.^2))*F;
            B=Best_vulture2_X-((Best_vulture2_X.*current_vulture_X)./(Best_vulture2_X-current_vulture_X.^2))*F;
            current_vulture_X=(A+B)/2;
        else
            current_vulture_X=random_vulture_X-abs(random_vulture_X-current_vulture_X)*F.*levyFlight(variables_no);
        end
    end
    % phase 2
    if  abs(F)>=0.5
        if rand<p3
            current_vulture_X=(abs((2*rand)*random_vulture_X-current_vulture_X))*(F+rand)-(random_vulture_X-current_vulture_X);
        else
            s1=random_vulture_X.* (rand()*current_vulture_X/(2*pi)).*cos(current_vulture_X);
            s2=random_vulture_X.* (rand()*current_vulture_X/(2*pi)).*sin(current_vulture_X);
            current_vulture_X=random_vulture_X-(s1+s2);
        end
    end
end

function [current_vulture_X] = exploration(current_vulture_X, random_vulture_X, F, p1, upper_bound, lower_bound)

    if rand<p1
        current_vulture_X=random_vulture_X-(abs((2*rand)*random_vulture_X-current_vulture_X))*F;
    else
        current_vulture_X=(random_vulture_X-(F)+rand()*((upper_bound-lower_bound)*rand+lower_bound));
    end
    
end

function [ o ]=levyFlight(d)
  
    beta=3/2;

    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;
    v=randn(1,d);
    step=u./abs(v).^(1/beta);

    o=step;

end

function [ X ] = boundaryCheck(X, lb, ub)

    for i=1:size(X,1)
            FU=X(i,:)>ub;
            FL=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
    end
end

