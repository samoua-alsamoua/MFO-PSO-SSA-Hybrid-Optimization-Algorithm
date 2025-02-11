function [bestSolution, bestFitness, iteration]=MFO_CS2(p, dimension, maxIteration)

settingAlg;

N=30;
Max_iteration=maxIteration;
lb = lbArray;
ub = ubArray;
dim=dimension;
pa=0.25;

%Initialize the positions of moths
Moth_pos=initialization(N,dim,ub,lb);
Moth_fitness=zeros(1,N);

Iteration=1;

% Main loop
%%MFO Stage
while Iteration<= Max_iteration/2
    
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
%     Best_flame_score=fitness_sorted(1);
%     Best_flame_pos=sorted_population(1,:);
      
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

%% CS Stage
nest = Moth_pos;
fitness = Moth_fitness;

[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,p);


while Iteration> Max_iteration/2
    
    if Iteration > Max_iteration
        break;
    end
        % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,lb,ub);   
     [~,~,nest,fitness]=get_best_nest(nest,new_nest,fitness,p);
    % Update the counter
      Iteration=Iteration+1; 
    % Discovery and randomization
      new_nest=empty_nests(nest,lb,ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,p);
    % Update the counter again
      Iteration=Iteration+1;
    % Find the best objective so far  
    if fnew<fmin
        fmin=fnew;
        bestnest=best;
    end
end
%%
bestSolution=bestnest;
bestFitness=fmin;
iteration=Iteration;
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


function nest=get_cuckoos(nest,best,Lb,Ub)
    n=size(nest,1);
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    for j=1:n
        s=nest(j,:);
        % This is a simple way of implementing Levy flights
        % For standard random walks, use step=1;
        % Levy flights by Mantegna's algorithm
        u=randn(size(s))*sigma;
        v=randn(size(s));
        step=u./abs(v).^(1/beta);
        % In the next equation, the difference factor (s-best) means that 
        % when the solution is the best solution, it remains unchanged.     
        stepsize=0.01*step.*(s-best);
        % Here the factor 0.01 comes from the fact that L/100 should the typical
        % step size of walks/flights where L is the typical lenghtscale; 
        % otherwise, Levy flights may become too aggresive/efficient, 
        % which makes new solutions (even) jump out side of the design domain 
        % (and thus wasting evaluations).
        % Now the actual random walks or flights
        s=s+stepsize.*randn(size(s));
        nest(j,:)=simplebounds(s,Lb,Ub);
    end
end

function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,p)
% Evaluating all new solutions
    for j=1:size(nest,1)
        fnew=testFunction(newnest(j,:)', p);
        if fnew<=fitness(j)
           fitness(j)=fnew;
           nest(j,:)=newnest(j,:);
        end
    end
    % Find the current best
    [fmin,K]=min(fitness) ;
    best=nest(K,:);
end

function new_nest=empty_nests(nest,Lb,Ub,pa)
    % A fraction of worse nests are discovered with a probability pa
    n=size(nest,1);
    % Discovered or not -- a status vector
    K=rand(size(nest))>pa;
    % In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
    % this cuckoo's egg is less likely to be discovered, thus the fitness should 
    % be related to the difference in solutions.  Therefore, it is a good idea 
    % to do a random walk in a biased way with some random step sizes.  
    % New solution by biased/selective random walks
    stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
    new_nest=nest+stepsize.*K;
    for j=1:size(new_nest,1)
        s=new_nest(j,:);
        new_nest(j,:)=simplebounds(s,Lb,Ub);  
    end
end

function s=simplebounds(s,Lb,Ub)
    % Apply the lower bound
    ns_tmp=s;
    I=ns_tmp<Lb;
    ns_tmp(I)=Lb(I);

    % Apply the upper bounds 
    J=ns_tmp>Ub;
    ns_tmp(J)=Ub(J);
    % Update this new move 
    s=ns_tmp;
end

