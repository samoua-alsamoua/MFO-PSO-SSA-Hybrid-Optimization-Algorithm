% % filepath='\Mac\Home\Desktop\MultilevelImageSegmentation(SamouaAlSamoua)\new\Yeni\image\test';
% 
% % Specify the folder where the images are located
% folder = '\\Mac\Home\Desktop\MultilevelImageSegmentation(SamouaAlSamoua)\new\Yeni\image\test';
% 
% % Get a list of all files in the folder with the specified extensions
% filePattern = fullfile(folder, '*.*');
% theFiles = dir(filePattern);
% 
% % Loop through each file and display its name and extension
% for k = 1 : length(theFiles)
%     baseFileName = theFiles(k).name;
%     [~, name, ext] = fileparts(baseFileName);
%     fprintf('%s%s\n', name, ext);
% end
% 

% Specify the folder where the images are located
folder = '\\Mac\Home\Desktop\MultilevelImageSegmentation(SamouaAlSamoua)\new\Yeni\image\test';

% Get a list of all files in the folder with the specified extensions
filePattern = fullfile(folder, '*.*');
theFiles = dir(filePattern);

% Loop through each file and display its name and extension with a comma
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    [~, name, ext] = fileparts(baseFileName);
    fprintf('''%s%s'',\n', name, ext);
end

{'1.jpg','101085.jpg','101087.jpg','102061.jpg','103070.jpg','105025.jpg','106024.jpg','108005.jpg','108070.jpg','108082.jpg','109053.jpg','119082.jpg','12084.jpg','123074.jpg','126007.jpg','130026.jpg','134035.jpg','14037.jpg','143090.jpg','145086.jpg','147091.jpg','148026.jpg','148089.jpg','156065.jpg','157055.jpg','159008.jpg','160068.jpg','16077.jpg','163085.jpg','167062.jpg','167083.jpg','170057.jpg','175032.jpg','175043.jpg','182053.jpg','189080.jpg','19021.jpg','196073.jpg','197017.jpg','208001.jpg','210088.jpg','21077.jpg','216081.jpg','219090.jpg','220075.jpg','223061.jpg','227092.jpg','229036.jpg','236037.jpg','24077.jpg','241004.jpg','241048.jpg','253027.jpg','253055.jpg','260058.jpg','271035.jpg','285079.jpg','291000.jpg','295087.jpg','296007.jpg','296059.jpg','299086.jpg','300091.jpg','302008.jpg','304034.jpg','304074.jpg','306005.jpg','33039.jpg','351093.jpg','361010.jpg','37073.jpg','376043.jpg','38082.jpg','38092.jpg','385039.jpg','41033.jpg','41069.jpg','42012.jpg','42049.jpg','43074.jpg','45096.jpg','54082.jpg','55073.jpg','58060.jpg','62096.jpg','65033.jpg','66053.jpg','69015.jpg','69020.jpg','69040.jpg','76053.jpg','78004.jpg','8023.jpg','85048.jpg','86000.jpg','86016.jpg','86068.jpg','87046.jpg','89072.jpg','97033.jpg'};




