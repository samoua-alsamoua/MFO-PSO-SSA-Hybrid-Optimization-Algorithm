function [segmented] = segment(image, levels)
    limites = [0 levels 255];
    tamanho = size(image);
    segmented(:,:) = image * 0;
    k = 1;
    for i = 1:tamanho(1, 1)
        for j = 1:tamanho(1, 2)
            while(k < size(limites, 2))
                if(image(i, j) >= limites(1, k) && image(i, j) <= limites(1, k + 1))
                    segmented(i,j,1)=limites(1,k);
                end
                k = k + 1;
            end
            k = 1;
        end
    end
end