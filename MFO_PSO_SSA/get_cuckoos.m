%% Get cuckoos
function nest=get_cuckoos(nest,best,Lb,Ub) 
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n
    s=nest(j,:);
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.001*step.*(s-best);
 
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);
end

% Application of simple constraints
function s=simplebounds(s,lb,ub)

Flag4ub=s>ub;
Flag4lb=s<lb;
s=s.*(~(Flag4ub+Flag4lb))+ub.*Flag4ub+lb.*Flag4lb;
  % Apply the lower bound
%   ns_tmp=s;
%   I=ns_tmp<Lb;
%   ns_tmp(I)=Lb(I);
%   
%   % Apply the upper bounds 
%   J=ns_tmp>Ub;
%   ns_tmp(J)=Ub(J);
%   % Update this new move 
%   s=ns_tmp;