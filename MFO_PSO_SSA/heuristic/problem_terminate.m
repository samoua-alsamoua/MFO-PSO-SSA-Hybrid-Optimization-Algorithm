function [SearchAgents_no, lb, ub] = problem_terminate(dimension)
    % Parameter settings:
    
    % SearchAgents_no (population) size
    SearchAgents_no = 50;

    % problem lower band 
    lowerBand = -100;
    n = dimension;
    lb = ones(1, n) * lowerBand;
    % problem upper band
    upperBand = 100;
    ub = ones(1, n) * upperBand;
end