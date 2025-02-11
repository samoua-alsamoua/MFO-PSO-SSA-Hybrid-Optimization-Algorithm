function [bestSolution, bestFitness, iteration] = MFO_PSO_SSA(p, dimension, maxIteration)

    settingAlg;

    % Parameters for MFO
    N = 30; % Number of moths
    Max_iteration = maxIteration;
    lb = lbArray;
    ub = ubArray;
    dim = dimension;

    % Initialize the positions of moths
    Moth_pos = initialization(N, dim, ub, lb);
    Moth_fitness = zeros(1, N);

    Iteration = 1;

    % MFO Stage (First 1/3 of iterations)
    while Iteration <= Max_iteration / 3

        % Number of flames (decreases over time)
        Flame_no = round(N - Iteration * ((N - 1) / Max_iteration));

        % Evaluate fitness of moths
        for i = 1:size(Moth_pos, 1)
            % Ensure moths stay within bounds
            Flag4ub = Moth_pos(i, :) > ub;
            Flag4lb = Moth_pos(i, :) < lb;
            Moth_pos(i, :) = (Moth_pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            % Calculate fitness
            Moth_fitness(1, i) = testFunction(Moth_pos(i, :)', p);
        end

        % Sort moths and update flames
        if Iteration == 1
            [fitness_sorted, I] = sort(Moth_fitness);
            sorted_population = Moth_pos(I, :);
            best_flames = sorted_population;
            best_flame_fitness = fitness_sorted;
        else
            double_population = [previous_population; best_flames];
            double_fitness = [previous_fitness, best_flame_fitness];
            [double_fitness_sorted, I] = sort(double_fitness);
            double_sorted_population = double_population(I, :);
            fitness_sorted = double_fitness_sorted(1:N);
            sorted_population = double_sorted_population(1:N, :);
            best_flames = sorted_population;
            best_flame_fitness = fitness_sorted;
        end

        % Store previous population and fitness
        previous_population = Moth_pos;
        previous_fitness = Moth_fitness;

        % Update moth positions
        a = -1 + Iteration * ((-1) / Max_iteration); % Linearly decreases from -1 to -2
        for i = 1:size(Moth_pos, 1)
            for j = 1:size(Moth_pos, 2)
                if i <= Flame_no
                    % Update position with respect to corresponding flame
                    distance_to_flame = abs(sorted_population(i, j) - Moth_pos(i, j));
                    b = 1;
                    t = (a - 1) * rand + 1;
                    Moth_pos(i, j) = distance_to_flame * exp(b .* t) .* cos(t .* 2 * pi) + sorted_population(i, j);
                else
                    % Update position with respect to one flame
                    distance_to_flame = abs(sorted_population(Flame_no, j) - Moth_pos(i, j));
                    b = 1;
                    t = (a - 1) * rand + 1;
                    Moth_pos(i, j) = distance_to_flame * exp(b .* t) .* cos(t .* 2 * pi) + sorted_population(Flame_no, j);
                end
            end
        end

        Iteration = Iteration + 1;
    end

    % PSO Stage (Second 1/3 of iterations)
    % Initialize PSO parameters
    nPop = N; % Population size (same as MFO)
    w = 0.9; % Inertia weight
    wdamp = 0.5 / Max_iteration; % Inertia weight damping ratio
    c1 = 2.0; % Personal learning coefficient
    c2 = 2.0; % Global learning coefficient
    VelMax = 0.1 * (ub - lb); % Velocity limits
    VelMin = -VelMax;

    % Initialize particles
    empty_particle.Position = [];
    empty_particle.Cost = [];
    empty_particle.Velocity = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];
    particle = repmat(empty_particle, nPop, 1);

    % Initialize global best
    GlobalBest.Cost = inf;

    % Initialize particles with MFO results
    for i = 1:nPop
        particle(i).Position = Moth_pos(i, :);
        particle(i).Velocity = zeros(1, dim);
        particle(i).Cost = Moth_fitness(i);
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update global best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end
    end

    % PSO main loop
    while Iteration <= (2 * Max_iteration) / 3
        for i = 1:nPop
            % Update velocity
            particle(i).Velocity = w * particle(i).Velocity ...
                + c1 * rand(1, dim) .* (particle(i).Best.Position - particle(i).Position) ...
                + c2 * rand(1, dim) .* (GlobalBest.Position - particle(i).Position);

            % Apply velocity limits
            particle(i).Velocity = max(particle(i).Velocity, VelMin);
            particle(i).Velocity = min(particle(i).Velocity, VelMax);

            % Update position
            particle(i).Position = particle(i).Position + particle(i).Velocity;

            % Apply position limits
            particle(i).Position = max(particle(i).Position, lb);
            particle(i).Position = min(particle(i).Position, ub);

            % Evaluate fitness
            particle(i).Cost = testFunction(particle(i).Position', p);

            % Update personal best
            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update global best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end
            end
        end

        % Update inertia weight
        w = w - wdamp;

        Iteration = Iteration + 1;
    end

    % SSA Stage
    % Initialize SSA parameters
    SalpPositions = zeros(N, dim); % Initialize salps with PSO results
    SalpFitness = zeros(1, N); % Initialize salp fitness with PSO results

    % Copy PSO results to SSA
    for i = 1:N
        SalpPositions(i, :) = particle(i).Position; % Copy positions
        SalpFitness(1, i) = particle(i).Cost; % Copy fitness values
    end

    % Sort salps and find the best solution
    [sorted_salps_fitness, sorted_indexes] = sort(SalpFitness);
    Sorted_salps = SalpPositions(sorted_indexes, :);
    FoodPosition = Sorted_salps(1, :); % Best salp (food source)
    FoodFitness = sorted_salps_fitness(1); % Best fitness

    % SSA main loop
    while Iteration <= Max_iteration
        c1 = 2 * exp(-(4 * Iteration / Max_iteration)^2); % Eq. (3.2) in the SSA paper

        for i = 1:size(SalpPositions, 1)
            SalpPositions = SalpPositions';

            if i <= N / 2
                % Update leader salps (Eq. 3.1 in the SSA paper)
                for j = 1:dim
                    c2 = rand();
                    c3 = rand();
                    if c3 < 0.5
                        SalpPositions(j, i) = FoodPosition(j) + c1 * ((ub(j) - lb(j)) * c2 + lb(j));
                    else
                        SalpPositions(j, i) = FoodPosition(j) - c1 * ((ub(j) - lb(j)) * c2 + lb(j));
                    end
                end
            elseif i > N / 2 && i < N + 1
                % Update follower salps (Eq. 3.4 in the SSA paper)
                point1 = SalpPositions(:, i - 1);
                point2 = SalpPositions(:, i);
                SalpPositions(:, i) = (point2 + point1) / 2;
            end

            SalpPositions = SalpPositions';
        end

        % Ensure salps stay within bounds and evaluate fitness
        for i = 1:size(SalpPositions, 1)
            Tp = SalpPositions(i, :) > ub;
            Tm = SalpPositions(i, :) < lb;
            SalpPositions(i, :) = (SalpPositions(i, :) .* (~(Tp + Tm))) + ub .* Tp + lb .* Tm;

            SalpFitness(1, i) = testFunction(SalpPositions(i, :)', p);

            % Update the food source (best solution)
            if SalpFitness(1, i) < FoodFitness
                FoodPosition = SalpPositions(i, :);
                FoodFitness = SalpFitness(1, i);
            end
        end

        Iteration = Iteration + 1;
    end

    % Return results
    bestSolution = FoodPosition;
    bestFitness = FoodFitness;
    iteration = Iteration;
end

function X=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end