function [bestSolution, bestFitness, iter]=LFD(pp, dimension, maxIteration)

	settingAlg;
    N = 50;
    lb = lbArray;
    ub = ubArray;
    threshold=2;
    dim = dimension;
    Positions=Initialization(N,dim, ub,lb);

    PositionsFitness = zeros(1,N);
    Positions_temp=Positions;

    for i=1:size(Positions,1)
        PositionsFitness(1,i) = testFunction(Positions(i,:)', pp);
    end
    [sorted_fitness,sorted_indexes]=sort(PositionsFitness);
    for newindex=1:N
        Sorted_Positions(newindex,:)=Positions(sorted_indexes(newindex),:);
    end
    TargetPosition=Sorted_Positions(1,:);
    TargetFitness=sorted_fitness(1);
    vec_flag=[1,-1];
    NN=[0,1];
    % Main loop
    l=1;
    while l<=maxIteration
        [m,ll]=sort(NN);
        for i=1:N
            S_i=zeros(1,dim);
            NeighborN=0;
            for j=1:N
                flag_index = floor(2*rand()+1);
                var_flag=vec_flag(flag_index);
                if i~=j
                    dis=Distance(Positions(i,:),Positions(j,:));
                    if (dis<threshold)
                        NeighborN=NeighborN+1;
                        D=(PositionsFitness(j)/(PositionsFitness(i)+eps));
                        D(NeighborN)=((.9*(D-min(D)))./(max(D(:))-min(D)+eps))+.1;
                        if l==2
                            rand_leader_index = floor(N*rand()+1);
                            X_rand = Positions(rand_leader_index, :);                        
                        else
    %                more opportunities for discovering unvisited pattern solutions
                            R=rand();
                            CSV=.5;
                            if R<CSV
                                rand_leader_index = floor(2*rand()+1);
                                X_rand = Positions(ll(rand_leader_index), :);
                                %Positions_temp(j,:)=Positions(j,:)+var_flag*.005*rand*(X_rand-Positions(j,:));
                                Positions_temp(j,:)=LevyFlights(Positions(j,:),X_rand,lb,ub);
                            else
                                Positions_temp(j,:)=lb(1)+rand(1,dim)*(ub(1)-lb(1));
                            end                       
                        end
                        pos_temp_nei{NeighborN}=Positions(j,:); 
                    end
                end
            end
            for p=1:NeighborN
                s_ij=var_flag*D(NeighborN).*(pos_temp_nei{p})/NeighborN;
                S_i=S_i+s_ij;
            end    
            S_i_total= S_i;
            rand_leader_index = floor(N*rand()+1);
            X_rand = Positions(rand_leader_index, :);    
            X_new = TargetPosition+10*S_i_total+rand*.00005*((TargetPosition+.005*X_rand)/2-Positions(i,:));
            X_new=LevyFlights(X_new,TargetPosition,lb,ub);
            Positions_temp(i,:)=X_new;
            NN(i)=NeighborN;
        end
        Positions=Positions_temp;
        PositionsFitness =  testFunction(Positions', pp);
        l = l + N;
        [xminn,x_pos_min]=min(PositionsFitness);
        if xminn<TargetFitness
            TargetPosition=Positions(x_pos_min,:);
            TargetFitness=xminn;
        end   
    end
 
    bestSolution = TargetPosition;
    bestFitness = TargetFitness;
    iter = l;

end

    function CP=LevyFlights(CP,DP,Lb,Ub)
        % Levy flights
        n=size(CP,1);
        % Levy exponent and coefficient
        % For details, see equation (2.21), Page 16 (chapter 2) of the book
        % X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
        beta=3/2;
        sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);     
        for j=1:n,
            s=CP(j,:);
            % This is a simple way of implementing Levy flights
            % For standard random walks, use step=1;
            %% Levy flights by Mantegna's algorithm
            u=randn(size(s))*sigma;
            v=randn(size(s));
            step=u./abs(v).^(1/beta);
            
            % In the next equation, the difference factor (s-best) means that
            % when the solution is the best solution, it remains unchanged.
            stepsize=0.01*step.*(s-DP);
            % Here the factor 0.01 comes from the fact that L/100 should the typical
            % step size of walks/flights where L is the typical lenghtscale;
            % otherwise, Levy flights may become too aggresive/efficient,
            % which makes new solutions (even) jump out side of the design domain
            % (and thus wasting evaluations).
            % Now the actual random walks or flights
            s=s+stepsize.*randn(size(s));
            % Apply simple bounds/limits
            CP(j,:)=simplebounds(s,Lb,Ub);
        end
    end

% Application of simple constraints
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
    
   
    function [X]=Initialization(N,dim,up,down)
        X = zeros(N,dim);
        for i=1:dim
            high=up(i);low=down(i);
            X(:,i)=rand(1,N).*(high-low)+low;
        end
    end
    
    function d = Distance(a,b)
        d=sqrt((a(1)-b(1))^2+(a(2)-b(2))^2);
    end
