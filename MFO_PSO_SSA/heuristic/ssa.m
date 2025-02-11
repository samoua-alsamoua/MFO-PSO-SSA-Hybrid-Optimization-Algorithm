function [bestSolution, bestFitness, iteration]=ssa(p, dimension, maxIteration)

settingAlg;

N=30;
Max_iter=ceil(maxIteration/N)+1;
dim=dimension;
lb = lbArray;
ub = ubArray;

%Initialize the positions of salps
SalpPositions=initialization(N,dim,ub,lb);

SalpFitness=zeros(1,N);
Sorted_salps=zeros(N,dim);

%calculate the fitness of initial salps

for i=1:size(SalpPositions,1)
    SalpFitness(1,i)=testFunction(SalpPositions(i,:)', p);
end

[sorted_salps_fitness,sorted_indexes]=sort(SalpFitness);

for newindex=1:N
    Sorted_salps(newindex,:)=SalpPositions(sorted_indexes(newindex),:);
end

FoodPosition=Sorted_salps(1,:);
FoodFitness=sorted_salps_fitness(1);

%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while l<Max_iter+1
    
    c1 = 2*exp(-(4*l/Max_iter)^2); % Eq. (3.2) in the paper
    
    for i=1:size(SalpPositions,1)
        
        SalpPositions= SalpPositions';
        
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    SalpPositions(j,i)=FoodPosition(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    SalpPositions(j,i)=FoodPosition(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>N/2 && i<N+1
            point1=SalpPositions(:,i-1);
            point2=SalpPositions(:,i);
            
            SalpPositions(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        SalpPositions= SalpPositions';
    end
    
    for i=1:size(SalpPositions,1)
        Tp=SalpPositions(i,:)>ub;
        Tm=SalpPositions(i,:)<lb;
        SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+ub.*Tp+lb.*Tm;
        
        SalpFitness(1,i)=testFunction(SalpPositions(i,:)', p);
        
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=SalpPositions(i,:);
            FoodFitness=SalpFitness(1,i);
            
        end
    end
    l = l + 1;
end
bestFitness=FoodFitness;
bestSolution=FoodPosition;
iteration=(l-2)*N;
end

function Positions=initialization(SearchAgents_no,dim,ub,lb)
    Positions = zeros(SearchAgents_no, dim);
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
  
end