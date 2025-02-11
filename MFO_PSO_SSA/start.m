clear;
clc;
addpath('heuristic');
run = 25;
level = 2;
% algorithms = {'pso','abc','cs','de','sos','ssa', 'sca','mfo','hho','gwo','woa'};
% algorithms = {'mfo','cs','MFO_CS1','MFO_CS2'};
algorithms = {'MFO_PSO_SSA'};

% images = {'1.jpg','101085.jpg','101087.jpg','102061.jpg','103070.jpg','105025.jpg','106024.jpg','108005.jpg','108070.jpg','108082.jpg','109053.jpg','119082.jpg','12084.jpg','123074.jpg','126007.jpg','130026.jpg','134035.jpg','14037.jpg','143090.jpg','145086.jpg','147091.jpg','148026.jpg','148089.jpg','156065.jpg','157055.jpg','159008.jpg','160068.jpg','16077.jpg','163085.jpg','167062.jpg','167083.jpg','170057.jpg','175032.jpg','175043.jpg','182053.jpg','189080.jpg','19021.jpg','196073.jpg','197017.jpg','208001.jpg','210088.jpg','21077.jpg','216081.jpg','219090.jpg','220075.jpg','223061.jpg','227092.jpg','229036.jpg','236037.jpg','24077.jpg','241004.jpg','241048.jpg','253027.jpg','253055.jpg','260058.jpg','271035.jpg','285079.jpg','291000.jpg','295087.jpg','296007.jpg','296059.jpg','299086.jpg','300091.jpg','302008.jpg','304034.jpg','304074.jpg','306005.jpg','33039.jpg','351093.jpg','361010.jpg','37073.jpg','376043.jpg','38082.jpg','38092.jpg','385039.jpg','41033.jpg','41069.jpg','42012.jpg','42049.jpg','43074.jpg','45096.jpg','54082.jpg','55073.jpg','58060.jpg','62096.jpg','65033.jpg','66053.jpg','69015.jpg','69020.jpg','69040.jpg','76053.jpg','78004.jpg','8023.jpg','85048.jpg','86000.jpg','86016.jpg','86068.jpg','87046.jpg','89072.jpg','97033.jpg'};
images = {'1.jpg'};
maxIteration = 500 * level; 

functionsNumber = length(images); 
out = zeros(length(algorithms), functionsNumber, run, 3, 4);
outSolution = zeros(length(algorithms), functionsNumber, 3, level);

for i = 1 : length(algorithms)
    disp(algorithms(i));
    algorithm = str2func(algorithms{i});
    for j = 1 : functionsNumber
        disp(j);
        image = imread(strcat('image/', images{j}));
        [~, ~, d] = size(image);
        if d == 1
            p = getGrayP(image);
            minFitness = inf;
            for k = 1 : run
                [bestSolution, bestFitness, ~] = algorithm(p, level, maxIteration);
                [m, ~] = size(bestSolution);
                if(m > 1) 
                    bestSolution = bestSolution';
                end
                if minFitness > bestFitness
                    minFitness = bestFitness;
                    outSolution(i, j, 1, :) = fix(sort(bestSolution));
                end
                bestFitness = 1 / bestFitness;
                segmentedImage = segment(image, fix(sort(bestSolution)));
                psnr = getPSNR(image, segmentedImage);
                ssim = getMSSIM(image, segmentedImage);
                fsim = getFSIM(image, segmentedImage);
                out(i, j, k, 1, :) = [bestFitness, psnr, ssim, fsim];
            end
        else
            for n = 1 : 3
                p = getColorP(image, n);
                minFitness = inf;
                for k = 1 : run
                    [bestSolution, bestFitness, ~] = algorithm(p, level, maxIteration);
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
                    out(i, j, k, n, :) = [bestFitness, psnr, ssim, fsim];
                end
            end
        end
    end
  
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(outSolution(i, :, 1, :)), 'Red-Solution');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 1)), 'Red-Fitness');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 2)), 'Red-PSNR');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 3)), 'Red-SSIM');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 1, 4)), 'Red-FSIM');

   
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(outSolution(i, :, 2, :)), 'Green-Solution');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 1)), 'Green-Fitness');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 2)), 'Green-PSNR');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 3)), 'Green-SSIM');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 2, 4)), 'Green-FSIM');


    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(outSolution(i, :, 3, :)), 'Blue-Solution');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 1)), 'Blue-Fitness');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 2)), 'Blue-PSNR');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 3)), 'Blue-SSIM');
    xlswrite(strcat(func2str(algorithm), '-d=', num2str(level), '.xlsx'), squeeze(out(i, :, :, 3, 4)), 'Blue-FSIM');
    
    eD = strcat(func2str(algorithm), '-Bitti :)');
    disp(eD);
end