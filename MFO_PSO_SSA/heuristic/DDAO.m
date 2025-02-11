function [bestSolution, bestFitness, iter]=DDAO(p, dimension, MaxIt)

settingAlg;
% Nvar = 10;              % Number of Variables
VarLength = [1 dimension];   % solution vector size
lb = lbArray;
ub = ubArray;
L_limit = lb;          % Lower limit of solution vector
U_limit = ub;           % Upper limit of solution vector

%% DDAO Parameters
Maxj=2 * MaxIt;    % Maximum Number of Sub-iterations
T0=2000;       % Initial Temp.
alpha=0.995;     % Temp. Reduction Rate
Npop=3;        % Population Size
%% Initialization
empty_template.Phase=[];
empty_template.Cost=[];
pop=repmat(empty_template,Npop,1);
% Initialize Best Solution
BestSol.Cost=inf;

%%initialization(N,dim,up,down)
% Initialize Population
for i=1:Npop    
    % Initialize Position
    pop(i).Phase= unifrnd(L_limit,U_limit,VarLength);
    % Evaluation
%             fitness(i) = testFunction(X(i, :)', p); % store fitness at iteration t        
    pop(i).Cost=testFunction(pop(i).Phase', p);    
    % Update Best Solution
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end    
end
% Vector to Hold Best Costs
BestCost=zeros(MaxIt,1);
% Intialize Temp.
T = T0;

t=0; % Loop counter
%% main loop
while t < MaxIt
    newpop = repmat(empty_template,Maxj,1);
    
        for j=1:Maxj        
        % Create and Evaluate New Solutions        
        newpop(j).Phase = unifrnd(L_limit,U_limit,VarLength);
        % set the new solution within the search space
        newpop(j).Phase = max(newpop(j).Phase, L_limit);
        newpop(j).Phase = min(newpop(j).Phase, U_limit);
        % Evaluate new solution
        newpop(j).Cost=testFunction(newpop(j).Phase', p);
        t = t + 1;
        end        
        % Sort Neighbors
        [~, SortOrder]=sort([newpop.Cost]);
        newpop=newpop(SortOrder);
        bnew = newpop(1);
        kk = randi(Npop);
        bb = randi(Npop);
        % forging parameter
        if(rem(t,2)==1)
        Mnew.Phase = (pop(kk).Phase-pop(bb).Phase)+ bnew.Phase;
        elseif (rem(t,2)==0)
             Mnew.Phase = (pop(kk).Phase-pop(bb).Phase)+ bnew.Phase*rand;
        end
        % set the new solution within the search space
        Mnew.Phase = max(Mnew.Phase, L_limit);
        Mnew.Phase = min(Mnew.Phase, U_limit); 
        % Evaluate new solution
        Mnew.Cost=testFunction(Mnew.Phase', p);
        t = t + 1;
        for i=1:Npop            
            if Mnew.Cost <= pop(i).Cost
                pop(i)= Mnew;                
            else
                DELTA=(Mnew.Cost-pop(i).Cost);
                P=exp(-DELTA/T);
                if rand <= P
                    pop(end)= Mnew;                                      
                end            
            end  
            % Update Best Solution Ever Found
            if pop(i).Cost <= BestSol.Cost
                BestSol=pop(i);
            end  
        end
     
    % Update Temp.
    T=alpha*T;
end

bestSolution = BestSol.Phase;
bestFitness = BestSol.Cost;
iter = t;

end

% function [X]=initialization(N,dim,up,down)
% 
% if size(up,1)==1
%     X=rand(N,dim).*(up-down)+down;
% end
% if size(up,1)>1
%     for i=1:dim
%         high=up(i);low=down(i);
%         X(:,i)=rand(1,N).*(high-low)+low;
%     end
% end
% end