function [bestSolution, bestFitness, iteration]=mfo(p, dimension, maxIteration)

settingAlg;

N=30;
Max_iteration=ceil(maxIteration/N)+1;
lb = lbArray;
ub = ubArray;
dim=dimension;

%Initialize the positions of moths
Moth_pos=initialization(N,dim,ub,lb);
Moth_fitness=zeros(1,N);

Iteration=1;

% Main loop
while Iteration<Max_iteration+1
    
    % Number of flames Eq. (3.14) in the paper
    Flame_no=round(N-Iteration*((N-1)/Max_iteration));
    
    for i=1:size(Moth_pos,1)
        
        % Check if moths go out of the search spaceand bring it back
        Flag4ub=Moth_pos(i,:)>ub;
        Flag4lb=Moth_pos(i,:)<lb;
        Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        % Calculate the fitness of moths
        Moth_fitness(1,i)=testFunction(Moth_pos(i,:)', p);
        
    end
       
    if Iteration==1
        % Sort the first population of moths
        [fitness_sorted, I]=sort(Moth_fitness);
        sorted_population=Moth_pos(I,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    else
        
        % Sort the moths
        double_population=[previous_population;best_flames];
        double_fitness=[previous_fitness best_flame_fitness];
        
        [double_fitness_sorted, I]=sort(double_fitness);
        double_sorted_population=double_population(I,:);
        
        fitness_sorted=double_fitness_sorted(1:N);
        sorted_population=double_sorted_population(1:N,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    end
    
    % Update the position best flame obtained so far
    Best_flame_score=fitness_sorted(1);
    Best_flame_pos=sorted_population(1,:);
      
    previous_population=Moth_pos;
    previous_fitness=Moth_fitness;
    
    % a linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a=-1+Iteration*((-1)/Max_iteration);
    
    for i=1:size(Moth_pos,1)
        
        for j=1:size(Moth_pos,2)
            if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
                
                % D in Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(i,j);
            end
            
            if i>Flame_no % Upaate the position of the moth with respct to one flame
                
                % Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(Flame_no,j);
            end
            
        end
        
    end
    
    Iteration=Iteration+1; 
end
bestSolution=Best_flame_pos;
bestFitness=Best_flame_score;
iteration=(Iteration-2)*N;
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

