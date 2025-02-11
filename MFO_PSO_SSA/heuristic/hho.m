function [bestSolution, bestFitness, iter]=hho(p, dimension, maxIteration)

settingAlg;
T = maxIteration;
N = 30;
dim = dimension;
lb = lbArray;
ub = ubArray;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
Fitness = testFunction(X', p); 
t=0; % Loop counter

while t < T
   [~, minIndex] = min(Fitness);
    Rabbit_Location=X(minIndex,:);
    E1=2*(1-(t/T)); % factor to show the decreaing energy of rabbit
    % Update the location of Harris' hawks
    for i=1:size(X,1)
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if q<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
            end
            X(i,:) = Check_boundries(X(i,:), ub, lb);
            Fitness(i) = testFunction(X(i,:)', p);
            t = t + 1;
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            if r>=0.5
                if abs(Escaping_Energy)<0.5 % Hard besiege
                    X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
                end

                if abs(Escaping_Energy)>=0.5  % Soft besiege
                    Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                    X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                end
                    X(i,:) = Check_boundries(X(i,:), ub, lb);
                    Fitness(i) = testFunction(X(i,:)', p);
                    t = t + 1;
            end

            %% phase 2: performing team rapid dives (leapfrog movements)
            
            if r<0.5
                if abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                    Jump_strength=2*(1-rand());
                    X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                    X1 = Check_boundries(X1, ub, lb);
                    FitX1 = testFunction(X1', p);
                    t = t + 1;
                    if FitX1<Fitness(i) % improved move?
                        X(i,:)=X1;
                        Fitness(i) = FitX1;
                    else % hawks perform levy-based short rapid dives around the rabbit
                        X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                        X2 = Check_boundries(X2, ub, lb);
                        FitX2 = testFunction(X2', p);
                        t = t + 1;
                        if (FitX2<Fitness(i)) % improved move?
                            X(i,:)=X2;
                            Fitness(i) = FitX2;
                            
                        end
                    end
                end

                if abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                    % hawks try to decrease their average location with the rabbit
                    Jump_strength=2*(1-rand());
                    X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                    X1 = Check_boundries(X1, ub, lb);
                    FitX1 = testFunction(X1', p);
                    t = t + 1;
                    if FitX1<Fitness(i) % improved move?
                        X(i,:)=X1;
                        Fitness(i) = FitX1;
                    else % Perform levy-based short rapid dives around the rabbit
                        X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                        X2 = Check_boundries(X2, ub, lb);
                        FitX2 = testFunction(X2', p);
                        t = t + 1;
                        if (FitX2<Fitness(i)) % improved move?
                            X(i,:)=X2;
                            Fitness(i) = FitX2;
                        end
                    end
                end
            end
        end
    end
end

bestSolution = X(minIndex,:);
bestFitness = Fitness(minIndex);
iter = t;
end

% ___________________________________
function o=Levy(d)
    beta=1.5;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
    o=step;
end

function X = Check_boundries(X, ub, lb)
        FU = X>ub;
        FL = X<lb;
        X =(X.*(~(FU+FL)))+ub.*FU+lb.*FL;
end

function [X]=initialization(N,dim,up,down)

    if size(up,2)==1
        X=rand(N,dim).*(up-down)+down;
    end
    if size(up,2)>1
        for i=1:dim
            high=up(i);low=down(i);
            X(:,i)=rand(1,N).*(high-low)+low;
        end
    end
end

