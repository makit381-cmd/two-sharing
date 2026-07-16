function PTDF = PTDF_matrix()
    %PTDF功率转移因子（简化求解）
    branches = table2array(readtable("C:\Users\19712\Desktop\new line data.xlsx",'Sheet','sheet1','Range','A4:B125'));
    G = digraph(branches(:,1),branches(:,2));
    PTDF = zeros(123,122);
    for i = 1:123
        if findnode (G,i) > 0
            path_nodes = shortestpath(G,1,i);
        else
            fprintf('error');
        end
        for j = 1:(length(path_nodes)-1)
            u = path_nodes(j);
            v = path_nodes(j+1);
            [is_found,line_idx] = ismember([u,v],branches,'rows');
            if ~is_found%如果没找到，翻一下from to
                [is_found, line_idx] = ismember([v, u], branches, 'rows');
            end
            if is_found
                PTDF(i,line_idx) = -1;
            end
        end
    end
end