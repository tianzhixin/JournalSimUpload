%% Please run the files below.

%%%%%% GraphTreeGeneration.m
%%%%%% MFD_Generation_ArbitraryNoRegions_n_max.m
%%%%%% Asys_ori_Creation.m

%% Main idea of this code:
% Calculate transfer flow at the critical point.

%% variables
% TCF is trip completion flow
% i=1:55, 55 is the number of regions.
% A is the adjacency matrix, with zero diagonal elements.
% nbi is the neighbor set of i, where nbi=find(A(i,:)==1);

% Q_sum_ni is the sum of critical TCF of the neighboring regions of i.
% Q_tf is a matrix of the transfer flow at the critical/selected point.
%      with Q_tf(i,j) is the transfer flow from i to j at the selected point.
% d_slt is a vector of the disturbance at the selected point.

% beta_slt is the beta_ij matrix at the selected point.
% beta_slt_all is the beta matrix at the selected point, with
% beta_slt_all(i,i) is nonzero.

% Asys is the system matrix.
% Bsys = B*K, where K is the distributed control law.
% k_fb is a 55 by 55 cell of state feedback gains, with each cell having two coefficients.
%       The first coefficient k_fb{i,j}(1) is the state feedback gain of
%       region j, which is the dominant region.
%       The second coefficient k_fb{i,j}(2) is the state feedback gain of
%       region i, which is the non-dominant region.


%% Calculate Q_tf (transfer flow at the critical/selected point)---%%

% The principle of Q_tf is:
% The outflow (exclude the network pheriphery flow) of a region cannot be greater than 60% of its TCF.
% The inflow (exclude the network pheriphery flow) of a region cannot be greater than 60% of its TCF.

function [Asys_ori,Bsys_ori,L]=Asys_ori_Creation(T,stp4,no_region,A_ori,A_adj_ori,G_dot_slt,G_slt,n_max,n_slt)
Q_tf=zeros(no_region);
Q_tf_tree=zeros(no_region);
rand('seed',stp4);
RandArray = rand(10000,1);
RandIndex = 1;
for i=1:no_region
    nbi=find(A_ori(i,:)==1);  % A_ori is the adjacency matrix of the network graph (0 diagonal) 
    Q_sum_ni=sum(G_slt(nbi));
    
    for j=1:length(nbi)
        nbj=find(A_ori(nbi(j),:)==1);
        Q_sum_nj=sum(G_slt(nbj));
        0.6*G_slt(i)*G_slt(nbi(j))/Q_sum_ni;  % The outflow (exclude the network pheriphery flow) of a region cannot be greater than 60% of its TCF.
        0.6*G_slt(nbi(j))*G_slt(i)/Q_sum_nj;  % The inflow (exclude the network pheriphery flow) of a region cannot be greater than 60% of its TCF.
        Q_tf(i,nbi(j))=(RandArray(RandIndex)*0.4+0.8)*min(0.6*G_slt(i)*G_slt(nbi(j))/Q_sum_ni,0.6*G_slt(nbi(j))*G_slt(i)/Q_sum_nj);
        RandIndex = RandIndex+1;
%         if A_adj_ori_cpc(i,nbi(j))==1
%             Q_tf(i,nbi(j)) = 0.15*Q_tf(i,nbi(j));    % A_adj_ori_cpc =  A_ori - A_adj_nd; Graph adjacency matrix - tree adjacency matrix
%         end
    end
end
% Q_tf_tree = Q_tf.*A_adj;


% dn=0 at the slt point. Thus the equilibrium equ is given by:

for i=1:no_region
    d_slt(i)=-0.2*G_slt(i)+sum(Q_tf(i,:))-sum(Q_tf(:,i)); % The inflow from the external region is 0.2*G_slt, because qei=0.4G_slt(i), beta_ei = 0.5;
%     d_slt_tree(i)=-0.2*G_slt(i)+sum(Q_tf_tree(i,:))-sum(Q_tf_tree(:,i)); % The inflow from the external region is 0.2*G_slt
end


% a must be less than 0, otherwise the outflow is greater than TCF.
% This is only for a check.
for i=1:no_region
    a(i)=sum(Q_tf(i,:))-d_slt(i)-G_slt(i);   % Why not a(i)=sum(Q_tf(i,:))-G_slt(i) ?
end

for i=1:no_region
    beta_slt(i,:)=Q_tf(i,:)./G_slt(i);
%     beta_slt_tree(i,:)=Q_tf_tree(i,:)./G_slt(i);
end


% beta_slt_all includes beta_ei, which controls the pheriphery flow.
beta_slt_all=beta_slt;
% beta_slt_all_tree=beta_slt_tree;

for i=1:no_region
    beta_slt_all(i,i)=0.5;
%     beta_slt_all_tree(i,i)=0.5;
end


%% Vectorize belt_slt_all and beta_slt_all_tree.

beta_slt_all_vector = zeros(size(find(beta_slt_all)));
% beta_slt_all_vector_tree = zeros(size(find(beta_slt_all_tree)));

rd_beta = 1;
for i = 1:no_region
    for j = 1:no_region
        if beta_slt_all(i,j)~=0
            beta_slt_all_vector(rd_beta) = beta_slt_all(i,j);
            rd_beta = rd_beta+1;
        end
    end
end

% rd_beta_tree = 1;
% for i = 1:no_region
%     for j = 1:no_region
%         if beta_slt_all_tree(i,j)~=0
%             beta_slt_all_vector_tree(rd_beta) = beta_slt_all_tree(i,j);
%             rd_beta_tree = rd_beta_tree+1;
%         end
%     end
% end
        
%% The system matrix calculation.
Asys_ori=zeros(no_region);
for i=1:no_region
    Asys_ori(i,i)=-G_dot_slt(i)*sum(beta_slt(i,:));
    nbi=find(A_ori(i,:)==1);
    for j=1:length(nbi)
        Asys_ori(i,nbi(j))=beta_slt(nbi(j),i)*G_dot_slt(nbi(j));
    end
end
L=zeros(no_region);
for i=1:no_region
    L(i,i)=1/(n_max(i)-n_slt(i));
end
AoriDT_nor = eye(no_region) + T*L*Asys_ori*inv(L);

%% Generate A_adj_tr_apr
% Asys_ori_abs = abs(Asys_ori);
% Asys_adj_ori = Asys_ori_abs + Asys_ori_abs';
% Asys_adj_ori = Asys_adj_ori - diag(diag(Asys_adj_ori));
% Asys_adj_ori = -Asys_adj_ori; % For minimum spanning tree
% 
% Gph_ori=graph(Asys_adj_ori);
% figure
% pic1 = plot(Gph_ori);
% [T_ori,pred] = minspantree(Gph_ori);
% highlight(pic1,T_ori)
% 
% figure;
% plot(T_ori);
% AdjT = adjacency(T_ori);
% A_adj = full(adjacency(T_ori)) + eye(no_region);  % A_adj_tr_apr is the adjacency matrix of the minimum spanning tree in Gph_ori
% Asys = Asys_ori.*A_adj;
% % A_adj = A_adj_tr_apr;
% % A_adj_nd = full(adjacency(T_ori));
% A_adj_ori_cpc = A_adj_ori - A_adj;


%% Tree network info
% Q_tf_tree = Q_tf.*A_adj;
% for i=1:no_region
% %     d_slt(i)=-0.2*G_slt(i)+sum(Q_tf(i,:))-sum(Q_tf(:,i)); % The inflow from the external region is 0.2*G_slt, because qei=0.4G_slt(i), beta_ei = 0.5;
%     d_slt_tree(i)=-0.2*G_slt(i)+sum(Q_tf_tree(i,:))-sum(Q_tf_tree(:,i)); % The inflow from the external region is 0.2*G_slt
% end
% for i=1:no_region
% %     beta_slt(i,:)=Q_tf(i,:)./G_slt(i);
%     beta_slt_tree(i,:)=Q_tf_tree(i,:)./G_slt(i);
% end
% 
% beta_slt_all_tree=beta_slt_tree;
% for i=1:no_region
% %     beta_slt_all(i,i)=0.5;
%     beta_slt_all_tree(i,i)=0.5;
% end
% beta_slt_all_vector_tree = zeros(size(find(beta_slt_all_tree)));
% rd_beta_tree = 1;
% for i = 1:no_region
%     for j = 1:no_region
%         if beta_slt_all_tree(i,j)~=0
%             beta_slt_all_vector_tree(rd_beta) = beta_slt_all_tree(i,j);
%             rd_beta_tree = rd_beta_tree+1;
%         end
%     end
% end

%% B matrix calculation
q_ei=0.4*G_slt;

% no_Inputs = length(find(A_adj));
% Bsys = zeros(no_region,no_Inputs);
% 
% index = 0;
% for i = 1:no_region
%     for j = 1:no_region
%         if A_adj(i,j)~=0
%             index = index + 1;
%             if i==j
%                 Bsys(i,index) = q_ei(i);
%             else
%                 Bsys(i,index) = -G_slt(i);
%                 Bsys(j,index) = G_slt(i);
%             end
%         end
%     end
% end

% Aadj_ori_triu = triu(A_adj_ori);
no_Inputs_ori = length(find(A_adj_ori));
Bsys_ori = zeros(no_region,no_Inputs_ori);

index = 0;
for i = 1:no_region
    for j = 1:no_region
        if A_adj_ori(i,j)~=0
            index = index + 1;
            if i==j
                Bsys_ori(i,index) = q_ei(i);
            else
                Bsys_ori(i,index) = -G_slt(i);
                Bsys_ori(j,index) = G_slt(i);
            end
        end
    end
end


% 
% clearvars -except no_region A_adj_ori  A_ori A_adj_ori_cpc A_adj...
%      G G_dot_slt G_slt n_max n_slt k_o k_n...
%      Asys Bsys Asys_ori beta_slt_all beta_slt_all_vector d_slt L...
%      Q_tf_tree d_slt_tree beta_slt_tree beta_slt_all_tree beta_slt_all_vector_tree

end

% clearvars -except  AdjT A_adj A_adj_nd no_region A_adj_ori A_adj_ori_cpc A_ori...
%      G G_dot_slt G_slt n_max n_slt k_o k_n...
%      Asys Bsys Asys_ori beta_slt_all beta_slt_all_vector d_slt L...
%      Q_tf_tree d_slt_tree beta_slt_tree beta_slt_all_tree beta_slt_all_vector_tree
 
% A_adj: Adjacency matrix of tree T_ori with 1 diagonal
% A_adj_ori: Adjacency matrix of G_ori, which is the perturbed T, with 1
%            diagonal
% A_adj_ori_cpc: A_adj_ori_cpc = ~A_adj_ori; The complement
% A_ori: Adjacency matrix of G_ori, which is the perturbed T, with 0
%            diagonal. A_ori = A_adj_ori - diag(diag(A_adj_ori));


% Asys: System matrix for the approximated tree network
% Asys_ori: System matrix for the original network
% Note that we should have A_adj_tr_apr(the approximation) = Asys(the original tree)
% beta_slt_all: beta hat matrix
% beta_slt_all_vector: vectorized beta hat
