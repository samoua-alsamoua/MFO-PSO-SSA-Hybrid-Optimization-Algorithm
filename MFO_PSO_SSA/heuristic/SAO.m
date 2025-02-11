%___________________________________________________________________%
%  Smell Agent Optimzation (SAO)              %

function [bestSolution, bestFitness, iter]=SAO(p, dimension, maxIteration)

settingAlg;
dim = dimension;
Agent_Pos=zeros(1,dim);
Agent_Fit=inf; %
lb = lbArray;
ub = ubArray;

N=50; % Number of Smell N (Search Agent)
olf=0.9;
K=0.6;
T=0.95;
M=0.9;
Step=0.02;

%Create the initial position of smell N
moles_Pos=initialization(N,dim,ub,lb);
fitness=testFunction(moles_Pos', p);
BestScore=inf;

l=0;% 

% Main loop
while l<maxIteration
    for i=1:size(moles_Pos,1)        
%Make Sure smell N remains in the search space.
        Clip_ub=moles_Pos(i,:)>ub;
        Clip_lb=moles_Pos(i,:)<lb;
        moles_Pos(i,:)=(moles_Pos(i,:).*(~(Clip_ub+Clip_lb)))+ub.*Clip_ub+lb.*Clip_lb;                      
        % Calculate objective function for each N
        fitness(i)=testFunction(moles_Pos(i, :)', p); 
        l=l+1;
        % Agent Fitness
        if fitness(i)<Agent_Fit 
            Agent_Fit=fitness(i); % Update Agent fitness
            Agent_Pos=moles_Pos(i,:);
        end
    end  
    % Update the Position of N
    for i=1:size(moles_Pos,1)
        for j=1:size(moles_Pos,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r3=rand();       
            r4=rand();
            r5=rand();
            Sniff_mole(i,j)=moles_Pos(i,j)+r1*sqrt(3*K*T/M); %Sniffing Mode
        end
        fitness(i)=testFunction(moles_Pos(i, :)', p); 
        l=l+1;
        [~,Index]=min(fitness(i));
        Agent_Pos=Sniff_mole(:,Index);
        [~,Indes]=max(fitness(i));
        Worst_Pos=Sniff_mole(:,Indes);
        if fitness(i)<BestScore
            BestScore=fitness(i);
            moles_Pos(i,:)=Sniff_mole(i,:);
        end
        %Trailing Mode       
        for j=1:size(moles_Pos,2)     
            Trail_mole(i,j)=moles_Pos(i,j)+r3*olf*(moles_Pos(i,j)-Agent_Pos(i,1))...
                -r4*olf*(moles_Pos(i,j)-Worst_Pos(i,1)); %Traili Mode
        end
        Trail_mole(i,:) = Check_boundries(Trail_mole(i,:), ub, lb);
        fitness(i)=testFunction(Trail_mole(i,:)', p); 
        l=l+1;
        if fitness(i)<BestScore
            BestScore=fitness(i);
            moles_Pos(i,:)=Trail_mole(i,:);
        end
         %Random Mode
        for j=1:size(moles_Pos,2)     
            Random_mole(i,j)=moles_Pos(i,j)+r5*Step; 
        end
        Random_mole(i,:) = Check_boundries(Random_mole(i,:), ub, lb);
        fitness(i)=testFunction(Random_mole(i,:)', p); 
        l=l+1;
        if fitness(i)<BestScore
            BestScore=fitness(i);
            moles_Pos(i,:)=Random_mole(i,:);
        end                  
    end
end
%     Agent_Fit=BestScore;
bestSolution = moles_Pos(1:dim);
bestFitness = Agent_Fit;
iter = l;
end

% This function initialize the first population of search agents
function Positions=initialization(nMol,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are are the same
if Boundary_no==1
    Positions=rand(nMol,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(nMol,1).*(ub_i-lb_i)+lb_i;
    end
end
end

function X = Check_boundries(X, ub, lb)
        FU = X>ub;
        FL = X<lb;
        X =(X.*(~(FU+FL)))+ub.*FU+lb.*FL;
end