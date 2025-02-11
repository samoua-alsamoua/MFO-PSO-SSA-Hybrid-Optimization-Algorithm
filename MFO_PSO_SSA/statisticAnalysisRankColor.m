clear;
clc;
algorithms = {'abc','aefa','aeo','agde', 'aso', 'boa','bsa','coa','cs','csa','de','dsa','efo','gwo','hho','is',...
    'lsa','lshade','mfla','mfo','mrfo','ms','pso','sca','sdo','se','sfs','sos','ssa','tlabc','wde','woa','yypo',...
    'fdbsos','fdbsfs','aoa','bes','bmo','cgo','choa','ddao','ecpo','ema','eo','fdbagde','fdbmrfo','fdbsdo'...
    'gbo','hfpso','hgso','igwo','lfd','lrfdbcoa','mao','mpso','msgo','mvo','pf','ppso','rso','sfo','spbo'...
    'sso','tfwo','tso'};
levels = [2 3 4 5 6 8 10 12 14 16];
algorithmNumber = length(algorithms);
imageCount = 15;
rank = 10;
source = 'kapur';
resultTable = cell(imageCount * length(levels), 8);
frequence = zeros(5, algorithmNumber);
color = zeros(1, imageCount * length(levels));
treshold = 0.001;

for i = 1 : length(levels)
    
    level = levels(i);

    file500 = strcat('500d\', source, '-result-d=', num2str(level), '.xlsx');
    file1000 = strcat('1000d\', source, '-result-d=', num2str(level), '.xlsx');
    file2000 = strcat('2000d\', source, '-result-d=', num2str(level), '.xlsx');
    file5000 = strcat('5000d\', source, '-result-d=', num2str(level), '.xlsx');

    data500 = xlsread(file500, 'Wilcoxon');
    data1000 = xlsread(file1000, 'Wilcoxon');
    data2000 = xlsread(file2000, 'Wilcoxon');
    data5000 = xlsread(file5000, 'Wilcoxon');
    
    fri500 = xlsread(file500, 'Friedman');
    fri1000 = xlsread(file1000, 'Friedman');
    fri2000 = xlsread(file2000, 'Friedman');
    fri5000 = xlsread(file5000, 'Friedman');
    
    friAll500 = fri500(:, 1:65) + fri500(:, 66:130) + fri500(:, 131:195);
    friAll1000 = fri1000(:, 1:65) + fri1000(:, 66:130) + fri1000(:, 131:195);
    friAll2000 = fri2000(:, 1:65) + fri2000(:, 66:130) + fri2000(:, 131:195);
    friAll5000 = fri5000(:, 1:65) + fri5000(:, 66:130) + fri5000(:, 131:195);
    
    disp(i);
    
    for j = 1 : imageCount
        [sorted500, indexes500] = sort(friAll500(j, :));
        [sorted1000, indexes1000] = sort(friAll1000(j, :));
        [sorted2000, indexes2000] = sort(friAll2000(j, :));
        [sorted5000, indexes5000] = sort(friAll5000(j, :));
        
        index500 = (j - 1) * algorithmNumber + indexes500(1:rank);
        index1000 = (j - 1) * algorithmNumber + indexes1000(1:rank);
        index2000 = (j - 1) * algorithmNumber + indexes2000(1:rank);
        index5000 = (j - 1) * algorithmNumber + indexes5000(1:rank);
            
        mean500 = mean([mean(data500(index500,2)) mean(data500(index500,8)) mean(data500(index500,14))]);
        std500 =  mean([mean(data500(index500,3)) mean(data500(index500,9)) mean(data500(index500,15))]);
        mean1000 = mean([mean(data1000(index1000,2)) mean(data1000(index1000,8)) mean(data1000(index1000,14))]);
        std1000 =  mean([mean(data1000(index1000,3)) mean(data1000(index1000,9)) mean(data1000(index1000,15))]);
        mean2000 = mean([mean(data2000(index2000,2)) mean(data2000(index2000,8)) mean(data2000(index2000,14))]);
        std2000 =  mean([mean(data2000(index2000,3)) mean(data2000(index2000,9)) mean(data2000(index2000,15))]);
        mean5000 = mean([mean(data5000(index5000,2)) mean(data5000(index5000,8)) mean(data5000(index5000,14))]);
        std5000 =  mean([mean(data5000(index5000,3)) mean(data5000(index5000,9)) mean(data5000(index5000,15))]);
        
        frequence(1, indexes500(1:rank)) = frequence(1, indexes500(1:rank)) + 1;
        frequence(2, indexes1000(1:rank)) = frequence(2, indexes1000(1:rank)) + 1;
        frequence(3, indexes2000(1:rank)) = frequence(3, indexes2000(1:rank)) + 1;
        frequence(4, indexes5000(1:rank)) = frequence(4, indexes5000(1:rank)) + 1;

        resultTable{(i - 1) * imageCount + j, 1} = strcat(num2str(mean500, '%10.6E'), ' (', num2str(std500, '%10.4E'), ')');
        resultTable{(i - 1) * imageCount + j, 2} = strcat(num2str(mean1000, '%10.6E'), ' (', num2str(std1000, '%10.4E'), ')');
        resultTable{(i - 1) * imageCount + j, 3} = strcat(num2str(mean2000, '%10.6E'), ' (', num2str(std2000, '%10.4E'), ')');
        resultTable{(i - 1) * imageCount + j, 4} = strcat(num2str(mean5000, '%10.6E'), ' (', num2str(std5000, '%10.4E'), ')');
        
        if( mean5000 - mean2000 < treshold)
            if( mean2000 - mean1000 < treshold)
                if( mean1000 - mean500 < treshold)
                    color((i - 1) * imageCount + j) = 1;
                else
                    color((i - 1) * imageCount + j) = 2;
                end
            else
                color((i - 1) * imageCount + j) = 3;
            end
        else
            color((i - 1) * imageCount + j) = 4;
        end
        
        resultTable{(i - 1) * imageCount + j, 5} = char(strjoin(algorithms(indexes500(1:rank))));
        resultTable{(i - 1) * imageCount + j, 6} = char(strjoin(algorithms(indexes1000(1:rank))));
        resultTable{(i - 1) * imageCount + j, 7} = char(strjoin(algorithms(indexes2000(1:rank))));
        resultTable{(i - 1) * imageCount + j, 8} = char(strjoin(algorithms(indexes5000(1:rank))));
        
    end
    
end

frequence(5, :) = sum(frequence(1:4, :));

filename = strcat('analiz-', source, '-rank-color=', num2str(rank), '.xlsx');
xlswrite(filename, resultTable);
pause(5);
sheetname = 'Sayfa1';
Excel = actxserver('excel.application');
WB = Excel.Workbooks.Open(fullfile('C:\Users\HP-Aras\Google Drive\1-Tolga-Sefa- Sezgisel Algoritmalar\2-Uygulamalar\segment workspace\', filename),0,false);
    
for i = 1 : (imageCount * length(levels))
    switch(color(i))
        case 1
            cell = strcat('A',num2str(i));
        case 2
            cell = strcat('B',num2str(i));
        case 3
            cell = strcat('C',num2str(i));
        case 4
            cell = strcat('D',num2str(i));
        otherwise
            cell = strcat('E',num2str(i));
    end
    WB.Worksheets.Item(1).Range(cell).Interior.ColorIndex = 4;
end
WB.Save();
WB.Close();
Excel.Quit();



