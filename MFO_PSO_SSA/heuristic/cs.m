function [bestSolution, bestFitness, iteration]=cs(p, dimension, maxIteration)

settingAlg;

pa=0.25;
n=25;
nd=dimension; 
Lb=lbArray; 
Ub=ubArray;

nest=zeros(n,nd);
for i=1:n
    nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
end
fitness=10^100*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,p);

N_iter=0;
while N_iter<maxIteration
    % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     [~,~,nest,fitness]=get_best_nest(nest,new_nest,fitness,p);
    % Update the counter
      N_iter=N_iter+n; 
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,p);
    % Update the counter again
      N_iter=N_iter+n;
    % Find the best objective so far  
    if fnew<fmin
        fmin=fnew;
        bestnest=best;
    end
end
bestSolution=bestnest;
bestFitness=fmin;
iteration=N_iter;
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

