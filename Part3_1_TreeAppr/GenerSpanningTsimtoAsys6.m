function [A_adj,AdjT,Asys,Bsys] = GenerSpanningTsimtoAsys3(NetworkNo,ptn,no_region,Asys_ori,G_slt)

rand('seed',NetworkNo);
% ptn = 2^(NetworkNo-4);
    

Asys_ori_pt = Asys_ori + Asys_ori.*(rand(no_region)*ptn);
Asys_ori_abs = abs(Asys_ori_pt);
Asys_adj_ori = Asys_ori_abs + Asys_ori_abs';
Asys_adj_ori = Asys_adj_ori - diag(diag(Asys_adj_ori));
Asys_adj_ori = -Asys_adj_ori; % For minimum spanning tree

Gph_ori=graph(Asys_adj_ori);
% figure
% pic1 = plot(Gph_ori);
[T_ori,pred] = minspantree(Gph_ori);
% highlight(pic1,T_ori)

% figure;
% plot(T_ori);
AdjT = adjacency(T_ori);
A_adj = full(adjacency(T_ori)) + eye(no_region);  % A_adj_tr_apr is the adjacency matrix of the minimum spanning tree in Gph_ori
Asys = Asys_ori.*A_adj;
% A_adj = A_adj_tr_apr;
% A_adj_nd = full(adjacency(T_ori));
% A_adj_ori_cpc = A_adj_ori - A_adj;

no_Inputs = length(find(A_adj));
Bsys = zeros(no_region,no_Inputs);
q_ei=0.4*G_slt;


index = 0;
for i = 1:no_region
    for j = 1:no_region
        if A_adj(i,j)~=0
            index = index + 1;
            if i==j
                Bsys(i,index) = q_ei(i);
            else
                Bsys(i,index) = -G_slt(i);
                Bsys(j,index) = G_slt(i);
            end
        end
    end
end

end