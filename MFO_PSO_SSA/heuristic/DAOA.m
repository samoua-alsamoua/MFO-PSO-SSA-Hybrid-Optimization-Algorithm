function [bestSolution, bestFitness, iter] = DAOA(p, dimension, maxIteration)

    settingAlg;
    Max_iter = maxIteration;
    dim = dimension;             % Number of Variable
    lb = lbArray;% Lower Bound
    ub = ubArray;% Upper Bound
    N = 5;                           %Number of Colors (Npop)
    Best_P = zeros(1, dim);           % bestSolution
    Best_FF = inf;                   % bestFitness
    Conv_curve = zeros(1, Max_iter);
    
    %% Initialize the positions of solution
    X = initialization(N, dim, ub, lb);
    Xnew = X;
    Ffun = zeros(1, size(X, 1));          % (fitness values)
    Ffun_new = zeros(1, size(Xnew, 1));   % (fitness values)

    l = 1;
    Mu = 0.001;
    alpha = 25;

    while l <= Max_iter                       % Main loop

        DAF = (Max_iter + 1 / l) ^ (alpha);    % DAF2
        DCS = 0.99 * ((1 - (l / Max_iter) ^ (0.5)));     % DCS

        % Update the Position of solutions
        for i = 1:size(X, 1)
            Ffun(i) = testFunction(X(i, :)', p); % Calculate the fitness values of solutions
            l = l + 1; 
            if Ffun(i) < Best_FF
                Best_FF = Ffun(i);
                Best_P = X(i, :);
            end

            for j = 1:size(X, 2)
                r1 = rand();
                if (size(lb, 2) == 1)
                    if r1 < DAF
                        r2 = rand();
                        if r2 > 0.5
                            Xnew(i, j) = (Best_P(1, j) / (DCS + eps) * ((ub - lb) * Mu + lb));
                        else
                            Xnew(i, j) = (Best_P(1, j) * DCS * ((ub - lb) * Mu + lb));
                        end
                    else
                        r3 = rand();
                        if r3 > 0.5
                            Xnew(i, j) = (Best_P(1, j) - DCS * ((ub - lb) * Mu + lb));
                        else
                            Xnew(i, j) = (Best_P(1, j) + DCS * ((ub - lb) * Mu + lb));
                        end
                    end
                end

                if (size(lb, 2) ~= 1)                          % if each of the ub and lb has more than one value
                    r1 = rand();
                    if r1 < DAF
                        r2 = rand();
                        if r2 > 0.5
                            Xnew(i, j) = ((Best_P(1, j) / (DCS + eps) * ((ub(j) - lb(j)) * Mu + lb(j))));
                        else
                            Xnew(i, j) = ((Best_P(1, j) * DCS * ((ub(j) - lb(j)) * Mu + lb(j))));
                        end
                    else
                        r3 = rand();
                        if r3 > 0.5
                            Xnew(i, j) = ((Best_P(1, j) - DCS * ((ub(j) - lb(j)) * Mu + lb(j))));
                        else
                            Xnew(i, j) = ((Best_P(1, j) + DCS * ((ub(j) - lb(j)) * Mu + lb(j))));
                        end
                    end
                end
            end

            Flag_ub = Xnew(i, :) > ub;                    % check if they exceed (up) the boundaries
            Flag_lb = Xnew(i, :) < lb;                    % check if they exceed (down) the boundaries
            Xnew(i, :) = (Xnew(i, :) .* (~(Flag_ub + Flag_lb))) + ub .* Flag_ub + lb .* Flag_lb;

            Ffun_new(i) = testFunction(Xnew(i, :)', p);          % calculate Fitness function
            l = l + 1;                                           % incremental iteration
            
            if Ffun_new(i) < Ffun(i)
                X(i, :) = Xnew(i, :);
                Ffun(i) = Ffun_new(i);
            end
            if Ffun(i) < Best_FF
                Best_FF = Ffun(i);
                Best_P = X(i, :);
            end
        end

        % Update the convergence curve
        Conv_curve(l) = Best_FF;

%         l = l + 1;  % incremental iteration
    end

    iter = l - 1;  % Correcting the iteration count
    bestSolution = Best_P;
    bestFitness = Best_FF;
end

function X = initialization(N, Dim, UB, LB)
    B_no = size(UB, 2); % numnber of boundaries
    X = zeros(N, Dim);

    % If each variable has a different lb and ub
    if B_no > 1
        for i = 1:Dim
            Ub_i = UB(i);
            Lb_i = LB(i);
            X(:, i) = rand(N, 1) .* (Ub_i - Lb_i) + Lb_i;
        end
    end
end
