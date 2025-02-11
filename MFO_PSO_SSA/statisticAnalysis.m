clear;
clc;
algorithms = {'WChOA','FDB_WChimp_CASE15','FDB_WChimp_CASE28','FDB_WChimp_CASE30','FDB_WChimp_CASE40','FDB_WChimp_CASE41', 'FDB_WChimp_CASE13','FDB_WChimp_CASE49','FDB_WChimp_CASE57', 'DAOA', 'DDAO', 'LFD', 'HHO', 'NCHHO', 'SAO', 'AVOA'};

levels = [2 4 6];
algorithmNumber = length(algorithms);
imageCount = 5;
source = '2000/';
resultTable = cell(imageCount * length(levels), 3);

for i = 1 : length(levels)
    
    level = levels(i);

    file500 = strcat('-result-d=', num2str(level), '.xlsx');
%     file1000 = strcat('1000d\', source, '-result-d=', num2str(level), '.xlsx');
%     file2000 = strcat('2000d\', source, '-result-d=', num2str(level), '.xlsx');

    data500 = xlsread(file500, 'Wilcoxon');
%     data1000 = xlsread(file1000, 'Wilcoxon');
%     data2000 = xlsread(file2000, 'Wilcoxon');
    
    disp(i);
    
    for j = 1 : imageCount
        fi = (j - 1) * algorithmNumber + 1;
        li = algorithmNumber * j;
    
        mean500 = mean([mean(data500(fi:li,2)) mean(data500(fi:li,8)) mean(data500(fi:li,14))]);
        std500 =  mean([mean(data500(fi:li,3)) mean(data500(fi:li,9)) mean(data500(fi:li,15))]);
%         mean1000 = mean([mean(data1000(fi:li,2)) mean(data1000(fi:li,8)) mean(data1000(fi:li,14))]);
%         std1000 =  mean([mean(data1000(fi:li,3)) mean(data1000(fi:li,9)) mean(data1000(fi:li,15))]);
%         mean2000 = mean([mean(data2000(fi:li,2)) mean(data2000(fi:li,8)) mean(data2000(fi:li,14))]);
%         std2000 =  mean([mean(data2000(fi:li,3)) mean(data2000(fi:li,9)) mean(data2000(fi:li,15))]);

        resultTable{(i - 1) * imageCount + j, 1} = strcat(num2str(mean500), '(', num2str(std500), ')');
%         resultTable{(i - 1) * imageCount + j, 2} = strcat(num2str(mean1000), '(', num2str(std1000), ')');
%         resultTable{(i - 1) * imageCount + j, 3} = strcat(num2str(mean2000), '(', num2str(std2000), ')');
    end
    
end

xlswrite('StaticAnalysis.xlsx', resultTable);



