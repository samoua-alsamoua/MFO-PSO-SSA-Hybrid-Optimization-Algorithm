clear;
clc;
global count;
addpath('heuristicF');
run = 25;
levels = [2 3 4 5 6 8 10 12 14 16];
algorithms = {'abc','aefa','aeo','agde', 'aso', 'boa','bsa','coa','cs','csa','de','dsa','efo','gwo','hho','is',...
    'lsa','lshade','mfla','mfo','mrfo','ms','pso','sca','sdo','se','sfs','sos','ssa','tlabc','wde','woa','yypo',...
    'fdbsos','fdbsfs','aoa','bes','bmo','cgo','choa','ddao','ecpo','ema','eo','fdbagde','fdbmrfo','fdbsdo'...
    'gbo','hfpso','hgso','igwo','lfd','lrfdbcoa','mao','mpso','msgo','mvo','pf','ppso','rso','sfo','spbo'...
    'sso','tfwo','tso'};
images = {'renkli/1.tiff','renkli/2.tiff','renkli/3.tiff','renkli/4.jpg','renkli/5.jpg','renkli/6.jpg','renkli/7.jpg',...
    'renkli/8.jpg','renkli/9.jpg','renkli/10.jpg','renkli/11.jpg','renkli/12.jpg','renkli/13.jpg','renkli/14.jpg','renkli/15.jpg'};

data = xlsread('fitness-kapur');

for li = 1 : 10
    level = levels(li);
    functionsNumber = length(images); 
    out = zeros(length(algorithms), functionsNumber, run, 3, 4);
    outSolution = zeros(length(algorithms), functionsNumber, 3, level);
    maxIteration = 55 * level;
    for i = 1 : length(algorithms)
        disp(algorithms(i));
        algorithm = str2func(algorithms{i});
        for j = 1 : functionsNumber
            disp(j);
            image = imread(strcat('image/', images{j}));
            for n = 1 : 3
                p = getColorP(image, n);
                minFitness = inf;
                for k = 1 : run
                    tic;
                    count = 0;
                    [bestSolution, bestFitness, success] = algorithm(p, level, maxIteration, 1 / data((i-1)*10 + j));
                    time = toc;
                    fe = count;
                    [m, ~] = size(bestSolution);
                    if(m > 1) 
                        bestSolution = bestSolution';
                    end
                    if minFitness > bestFitness
                        minFitness = bestFitness;
                        outSolution(i, j, n, :) = fix(sort(bestSolution));
                    end
                    bestFitness = 1 / bestFitness;
                    sImage = image(:, :, n);
                    segmentedImage = segment(sImage, fix(sort(bestSolution)));
                    psnr = getPSNR(sImage, segmentedImage);
                    ssim = getMSSIM(sImage, segmentedImage);
                    fsim = getFSIM(sImage, segmentedImage);
                    out(i, j, k, n, :) = [bestFitness, fe, time, success];
                end
            end
        end
        
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(outSolution(i, :, 1, :)), 'Red-Solution');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 1)), 'Red-Fitness');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 2)), 'Red-FE');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 3)), 'Red-Time');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 4)), 'Red-Success');

        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(outSolution(i, :, 2, :)), 'Green-Solution');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 1)), 'Green-Fitness');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 2)), 'Green-FE');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 3)), 'Green-Time');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 4)), 'Green-Success');

        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(outSolution(i, :, 3, :)), 'Blue-Solution');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 1)), 'Blue-Fitness');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 2)), 'Blue-FE');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 3)), 'Blue-Time');
        xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 4)), 'Blue-Success');

        eD = strcat(func2str(algorithm), '-Bitti :)');
        disp(eD);
    end
end