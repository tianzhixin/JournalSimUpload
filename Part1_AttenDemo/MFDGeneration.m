%% Please run the files below.

%%%%%% GraphTreeGeneration.m
%%%%%% MFD_Generation_ArbitraryNoRegions_n_max.m
%%%%%% Asys_ori_Creation.m

function [G_dot_slt,G_slt,n_max,n_slt] = MFD_Generation_ArbitraryNoRegions_n_max(no_region,stp3)

%% 
a_pt = 0.8;
nx=zeros(no_region,500);

for i=1:no_region
    for j=1:500
        nx(i,j)=10*(j-1);
    end
end

%%%%%%%%%%%------------original MFD of Yokohama-----------------%%%%%%%%
% 
% syms n0
% G0=((1.4877e-7)*(n0^3)-(2.9815e-3)*(n0^2)+15.0912*n0)/3600; %original MFD

%%%%%%%%%%-----------New MFD for each region-----------------%%%%%%%%
% To produce different MFDs, we use Yokohama MFD as the basic model, and
% scale it with different proportions.
p_min=0.3;  %p_min is the minimun proportion of the original MFD
p_max=0.5;   %p_max is the maximum proportion of the original MFD

rand('seed',stp3);
RandArray = rand(no_region,2);

for i=1:no_region

    k_o(i)=p_min+RandArray(i,1)*(p_max-p_min); %k_o is the proportion of the trip completion flow of the original MFD
    k_n(i)=p_min+RandArray(i,2)*(p_max-p_min); %k_n is the proportion of the accumulation of the original MFD

    n_slt(i)=a_pt*k_n(i)*3400;  %3400 is the accumulation corresponding to trip completion flow of Yokohama MFD
    n_max(i)=k_n(i)*10000;  %10000 is the max accumulation of Yokohama MFD
    G_slt(i)=k_o(i)*(((1.4877e-7)*((n_slt(i)./k_n(i)).^3)-(2.9815e-3)*((n_slt(i)./k_n(i)).^2)+15.0912*(n_slt(i)./k_n(i)))./3600);

    n_slt_u(i)=a_pt*k_n(i)*3400+0.1;
    n_slt_l(i)=a_pt*k_n(i)*3400-0.1;
    G_slt_u(i)=k_o(i)*(((1.4877e-7)*((n_slt_u(i)./k_n(i)).^3)-(2.9815e-3)*((n_slt_u(i)./k_n(i)).^2)+15.0912*(n_slt_u(i)./k_n(i)))./3600);
    G_slt_l(i)=k_o(i)*(((1.4877e-7)*((n_slt_l(i)./k_n(i)).^3)-(2.9815e-3)*((n_slt_l(i)./k_n(i)).^2)+15.0912*(n_slt_l(i)./k_n(i)))./3600);
    G_dot_slt(i)=(G_slt_u(i)-G_slt_l(i))/(n_slt_u(i)-n_slt_l(i));    
    
    
%     n(i)=n0/k_n(i);
    n_max(i)=10000*k_n(i);
    index_n_0{i}=find(nx(i,:)>n_max(i));
    nx(i,find(nx(i,:)>n_max(i)))=0;
    
    G(i,:)=k_o(i)*(((1.4877e-7)*((nx(i,:)./k_n(i)).^3)...
        -(2.9815e-3)*((nx(i,:)./k_n(i)).^2)+15.0912*(nx(i,:)./k_n(i)))./3600);
end

%% ------------------Plot MFD------------------

% figure

% for i=1:no_region
%      plot(nx(i,1:(min(index_n_0{i})-1)),G(i,1:(min(index_n_0{i})-1)))
%      hold on
%      plot(n_slt(i),G_slt(i),'*')
% end

% clearvars -except  AdjT A_adj A_adj_nd no_region A_adj_ori A_adj_ori_cpc A_ori...
%      G G_dot_slt G_slt n_max n_slt k_o k_n
end
     
 
% A_adj: Adjacency matrix of tree T_ori with 1 diagonal
% A_adj_ori: Adjacency matrix of G_ori, which is the perturbed T, with 1
%            diagonal
% A_adj_ori_cpc: A_adj_ori_cpc = ~A_adj_ori; The complement
% A_ori: Adjacency matrix of G_ori, which is the perturbed T, with 0
%            diagonal. A_ori = A_adj_ori - diag(diag(A_adj_ori));