
%__________________________________________________________________ %
%                          Multi-Objective                          %
%                Multi-objective Artificial Vultures                %
%                    Optimization Algorithm (MOAVOA)                %



function [BestSolution, bestFitness, iter] = MOAVOA(p, dimension, maxIteration)




% if method==3
% 
%     TestProblem=sprintf('P%d',m);
% 
%     fobj= Ptest(TestProblem);
% 
%     xrange  = xboundaryP(TestProblem);
%     dim=max(size(xrange));
%   
%     lb=xrange(:,1)';
%     ub=xrange(:,2)';
% 
% end
%

% Repository Size
settingAlg;
Alpha=0.1;  % Grid Inflation Parameter
nGrid=10;   % Number of Grids per each Dimension
Beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;
N = 100;        %%%% N
Archive_size = N;
Max_iter = maxIteration;
dim = dimension;
obj_no = 5;
lb = lbArray;
ub = ubArray;
X=CreateEmptyParticle(N);
current_vulture_X=CreateEmptyParticle(N);
Best_vulture1_X=CreateEmptyParticle(N);
Best_vulture2_X=CreateEmptyParticle(N);
% Initialize best solution and fitness
bestSolution.Position = zeros(1, dimension);
bestSolution.Cost = inf;

for i=1:N
    current_vulture_X(i,:).Position=zeros(1,dim);
    Best_vulture1_X(i,:).Position=zeros(1,dim);
    Best_vulture2_X(i,:).Position=zeros(1,dim);
    current_vulture_X(i,:).Cost=zeros(1,obj_no);
    Best_vulture1_X(i,:).Cost=inf(1,obj_no);
    Best_vulture2_X(i,:).Cost=inf(1,obj_no);
    current_vulture_X(i,:).Velocity=0;
    X(i,:).Velocity=0;
    X(i,:).Position=zeros(1,dim);


    current_vulture_X(i,:).Best.Position=current_vulture_X(i,:).Position;
    X(i,:).Best.Position=X(i,:).Position;
    current_vulture_X(i,:).Best.Cost=X(i,:).Cost;
    X(i,:).Best.Cost=X(i,:).Cost;
    X(i,:).Position=initialization(N,dim,ub,lb);
    X(i,:).Cost=testFunction(X(i,:).Position', p);
end


X=DetermineDominations(X);
Archive=GetNonDominatedParticles(X);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,Alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end


%%  Controlling parameter
p1=0.6;
p2=0.4;
p3=0.6;
alpha=0.8;
betha=0.2;
Gamma=2.5;

l = 0; % Loop counter

%%Main loop

while l < Max_iter
    Leader=SelectLeader(Archive,Beta);

    a=unifrnd(-2,2,1,1)*((sin((pi/2)*(l/Max_iter))^Gamma)+cos((pi/2)*(l/Max_iter))-1);
    P1=(2*rand+1)*(1-(l/Max_iter))+a;

    % ubdate the location
    for i=1:N
        current_vulture_X(i).Position = X(i).Position;  % pick the current vulture back to the population
        F=P1*(2*rand()-1);

        random_vulture_X(i).Position=random_select(current_vulture_X(i).Position, Leader.Position,alpha,betha);

        if abs(F) >= 1 % Exploration:
            current_vulture_X(i).Position = exploration(current_vulture_X(i).Position,random_vulture_X(i).Position, F, p1, ub, lb);
        elseif abs(F) < 1 % Exploitation:
            current_vulture_X(i).Position= exploitation(current_vulture_X(i).Position,  random_vulture_X(i).Position,Leader.Position, random_vulture_X(i).Position, F, p2, p3, dim, ub, lb);
        end
% ubdate best solution and fitness
        if current_vulture_X(i).Cost < bestSolution.Cost
            bestSolution.Position = current_vulture_X(i).Position;
            bestSolution.Cost = current_vulture_X(i).Cost;
        end
        
        current_vulture_X(i).Position= boundaryCheck( current_vulture_X(i).Position, lb, ub);

        current_vulture_X(i).Cost=testFunction(current_vulture_X(i).Position', p);
        l = l + 1;
        
    end

    current_vulture_X(i).Position= boundaryCheck( current_vulture_X(i).Position, lb, ub);

    current_vulture_X(i).Cost=testFunction(current_vulture_X(i).Position', p);
    l = l + 1;

    current_vulture_X=DetermineDominations(current_vulture_X);
    non_dominated_current_vulture_X=GetNonDominatedParticles(current_vulture_X);

    Archive=[Archive
        non_dominated_current_vulture_X];

    Archive=DetermineDominations(Archive);
    Archive=GetNonDominatedParticles(Archive);


    for i=1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end

    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);

        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,Alpha);

    end

%     disp(['In iteration ' num2str(l) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
%     save results

    % Results

    Archive_F=GetCosts(Archive);
end
bestFitness = bestSolution.Cost;
BestSolution = bestSolution.Position;
iter = l;
end

%%%%%% other functions
function particle=CreateEmptyParticle(n)
    
    if nargin<1
        n=1;
    end

    empty_particle.Position=[];
    empty_particle.Velocity=[];
    empty_particle.Cost=[];
    empty_particle.Dominated=false;
    empty_particle.Best.Position=[];
    empty_particle.Best.Cost=[];
    empty_particle.GridIndex=[];
    empty_particle.GridSubIndex=[];
    
    particle=repmat(empty_particle,n,1);
    
end

function G=CreateHypercubes(costs,ngrid,alpha)

    nobj=size(costs,1);
    
    empty_grid.Lower=[];
    empty_grid.ubper=[];
    G=repmat(empty_grid,nobj,1);
    
    for j=1:nobj
        
        min_cj=min(costs(j,:));
        max_cj=max(costs(j,:));
        
        dcj=alpha*(max_cj-min_cj);
        
        min_cj=min_cj-dcj;
        max_cj=max_cj+dcj;
        
        gx=linspace(min_cj,max_cj,ngrid-1);
        
        G(j).Lower=[-inf gx];
        G(j).ubper=[gx inf];
        
    end

end

function [ X ] = boundaryCheck(X, lb, ub)

    for i=1:size(X,1)
            FU=X(i,:)>ub;
            FL=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
    end
end

function pop=DetermineDominations(pop)

    npop=numel(pop);
    
    for i=1:npop
        pop(i).Dominated=false;
        for j=1:i-1
            if ~pop(j).Dominated
                if Dominates(pop(i),pop(j))
                    pop(j).Dominated=true;
                elseif Dominates(pop(j),pop(i))
                    pop(i).Dominated=true;
                    break;
                end
            end
        end
    end

end

function dom=Dominates(x,y)

    if isstruct(x)
        x=x.Cost;
    end

    if isstruct(y)
        y=y.Cost;
    end
    
    dom=all(x<=y) && any(x<y);

end

function [current_vulture_X] = exploitation(current_vulture_X, Best_vulture1_X, Best_vulture2_X, ...
                                                                      random_vulture_X, F, p2, p3, variables_no, ubper_bound, lower_bound)

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
            current_vulture_X=(abs(2*rand)*random_vulture_X-current_vulture_X)*(F+rand)-(random_vulture_X-current_vulture_X);
        else
            s1=random_vulture_X.* (rand()*current_vulture_X/(2*pi)).*cos(current_vulture_X);
            s2=random_vulture_X.* (rand()*current_vulture_X/(2*pi)).*sin(current_vulture_X);
            current_vulture_X=random_vulture_X-(s1+s2);
        end
    end
end

function [current_vulture_X] = exploration(current_vulture_X, random_vulture_X, F, p1, ubper_bound, lower_bound)

    if rand<p1
        current_vulture_X=random_vulture_X-(abs(2*rand)*random_vulture_X-current_vulture_X)*F;
    else
        current_vulture_X=(random_vulture_X-(F)+rand()*((ubper_bound-lower_bound)*rand+lower_bound));
    end
    
end

function nd_pop=GetNonDominatedParticles(pop)

    ND=~[pop.Dominated];
    
    nd_pop=pop(ND);

end

function costs=GetCosts(pop)

    nobj=numel(pop(1).Cost);
    costs=reshape([pop.Cost],nobj,[]);

end

function [Index, SubIndex]=GetGridIndex(particle,G)

    c=particle.Cost;
    
    nobj=numel(c);
    ngrid=numel(G(1).ubper);
    
    str=['sub2ind(' mat2str(ones(1,nobj)*ngrid)];

    SubIndex=zeros(1,nobj);
    for j=1:nobj
        
        U=G(j).ubper;
        
        i=find(c(j)<U,1,'first');
        
        SubIndex(j)=i;
        
        str=[str ',' num2str(i)];
    end
    
    str=[str ');'];
    
    Index=eval(str);
    
end

function [occ_cell_index, occ_cell_member_count]=GetOccubiedCells(pop)

    GridIndices=[pop.GridIndex];
    
    occ_cell_index=unique(GridIndices);
    
    occ_cell_member_count=zeros(size(occ_cell_index));

    m=numel(occ_cell_index);
    for k=1:m
        occ_cell_member_count(k)=sum(GridIndices==occ_cell_index(k));
    end
    
end

function rep_h=SelectLeader(rep,beta)
    if nargin<2
        beta=1;
    end

    [occ_cell_index, occ_cell_member_count]=GetOccubiedCells(rep);
    
    p=occ_cell_member_count.^(-beta);
    p=p/sum(p);
    
    selected_cell_index=occ_cell_index(RouletteWheelSelection(p));
    
    GridIndices=[rep.GridIndex];
    
    selected_cell_members=find(GridIndices==selected_cell_index);
    
    n=numel(selected_cell_members);
    
    selected_memebr_index=randi([1 n]);
    
    h=selected_cell_members(selected_memebr_index);
    
    rep_h=rep(h);
end

function i=RouletteWheelSelection(p)

    r=rand;
    c=cumsum(p);
    i=find(r<=c,1,'first');

end


function [random_vulture_X]=random_select(Best_vulture1_X,Best_vulture2_X,alpha,betha)

    probabilities=[alpha, betha ];
    
    if (RouletteWheelSelection( probabilities ) == 1)
            random_vulture_X=Best_vulture1_X;
    else
            random_vulture_X=Best_vulture2_X;
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

function rep=DeleteFromRep(rep,EXTRA,gamma)

    if nargin<3
        gamma=1;
    end

    for k=1:EXTRA
        [occ_cell_index, occ_cell_member_count]=GetOccupiedCells(rep);

        p=occ_cell_member_count.^gamma;
        p=p/sum(p);

        selected_cell_index=occ_cell_index(RouletteWheelSelection(p));

        GridIndices=[rep.GridIndex];

        selected_cell_members=find(GridIndices==selected_cell_index);

        n=numel(selected_cell_members);

        selected_memebr_index=randi([1 n]);

        j=selected_cell_members(selected_memebr_index);
        
        rep=[rep(1:j-1); rep(j+1:end)];
    end
    
end

function [occ_cell_index, occ_cell_member_count]=GetOccupiedCells(pop)

    GridIndices=[pop.GridIndex];
    
    occ_cell_index=unique(GridIndices);
    
    occ_cell_member_count=zeros(size(occ_cell_index));

    m=numel(occ_cell_index);
    for k=1:m
        occ_cell_member_count(k)=sum(GridIndices==occ_cell_index(k));
    end
    
end

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