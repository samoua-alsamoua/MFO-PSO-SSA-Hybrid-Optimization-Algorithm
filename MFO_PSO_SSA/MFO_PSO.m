function [bestSolution, bestFitness, iteration] = MFO_PSO(p, dimension, maxIteration)

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

    % MFO Stage
    while Iteration <= Max_iteration / 2

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

    % PSO Stage
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

    % Initialize particles with MFO results (The best solutions from MFO are passed to the PSO stage)
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
    while Iteration <= Max_iteration
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
            Iteration = Iteration + 1;


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

    end

    % Return results
    bestSolution = GlobalBest.Position;
    bestFitness = GlobalBest.Cost;
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