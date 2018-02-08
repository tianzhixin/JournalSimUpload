%% This code is to generate the adjacency matrix for a tree and the original graph.
%
%% The requirements are
% 1) the original graph must be a planar graph, such that it is a
% reasonable representation for the region partitioning in perimeter
% control.
% 2) The maximum spanning tree is more intense in the center area than the
% periphery area, for the same reason as above.
%
%% Main idea:
% 1) Generate a generally connected network by using the bubble model published by Kaixiang Wang, see Zhihu.
% 2) Assign weightings to edges. Particularly, the closer to center, the
% lighter the bubble is. Generate graph G.
% 3) Find the minimum spanning tree T in G.
%
%% Detailed explanation
% From "BubbleSimulatorPublish.m", we have 
% "BubbleResult50_main.mat" and "BubbleResult50.mat", which include
% 1) position info of those bubbles (posx, posy), 
% 2) overlapDist (overlapDist = bubbles.sumRadius-relDist, where
%                           relDist=sqrt(relPosx.^2+relPosy.^2);
% 3) N=overlapDist>0
%% How to get BubbleResult50_main.mat and BubbleResult50.mat?
% Run bubbleSimulatorPublish.m. See Kaixiang Wang's zhihu account for
% details.

%% Please run the files below.

%%%%%% GraphTreeGeneration.m
%%%%%% MFD_Generation_ArbitraryNoRegions_n_max.m
%%%%%% Asys_ori_Creation.m

%% Code.

% clear

% close all
function [A_adj_ori,A_ori] = GraphTreeGeneration(numBubbles,posx,posy,overlapDist)
% circle(posx,posy,ones(numBubbles,1)*100)    %Draw circles at position (posx, posy) with radius = 100
% axis equal
% rand('seed',stp2)

% for i=1:numBubbles
%     for j=1:numBubbles
%         if N(i,j)                                           %N=overlapDist>0; 
%         hold on, plot([posx(i),posx(j)],[posy(i),posy(j)])          % If circle i and circle j overlap, then draw an edge connecting them.
%         end
%     end
% end
Dist2Center = (posx.^2+posy.^2)/10e7;

wt = overlapDist;
% wt = rand(numBubbles).*overlapDist;
wt(wt<=0) = 0;               %There is no edge (thus no weight) between those non-touching circles.
wt_mf = zeros(size(wt));

wt_mf = wt;
wt_mf = wt_mf+wt_mf';
% figure
G=graph(wt_mf);
% p = plot(G);
% [T,pred] = minspantree(G);
% highlight(p,T)
%% dd
% createfigure(T)


%% dd
% figure;

% plot(T);
AdjG = adjacency(G);
% AdjT = adjacency(T);
% A_adj_nd = full(AdjT);
% A_adj = A_adj_nd + eye(numBubbles);
% no_region = numBubbles;
A_ori = full(AdjG);
A_adj_ori = A_ori + eye(numBubbles);
% A_adj_ori_cpc = A_ori - A_adj_nd;

%% Plot the tree graph in a planar gragh way
% figure
% for i=1:no_region
%     for j=1:i
%         if A_adj_nd(i,j)                                           %N=overlapDist>0; 
%         hold on, plot([posx(i),posx(j)],[posy(i),posy(j)],'LineWidth',1.5,'color',[0 0.5 1])          % If circle i and circle j overlap, then draw an edge connecting them.
%         end
%     end
% end
% 
% plot(posx,posy,'.','Markersize',30,'color',[0 0 1])
% 
% for i=1:no_region
%     txt{i}=num2str(i); 
% end



% for i = 1:no_region
% %     if i==47
% %         text(posx(i)-180,posy(i)-80,txt{i},'FontSize',13)
% %     else
% %         text(posx(i)+40,posy(i)+40,txt{i},'FontSize',13)
% %     end
%         text(posx(i)+40,posy(i)+40,txt{i},'FontSize',13)
% end

end
% clearvars -except AdjG AdjT A_adj A_adj_nd no_region A_adj_ori A_adj_ori_cpc A_ori

% A_adj: Adjacency matrix of tree T_ori with 1 diagonal
% A_adj_ori: Adjacency matrix of G_ori, which is the perturbed T, with 1
%            diagonal
% A_adj_ori_cpc: A_adj_ori_cpc = ~A_adj_ori; The complement
% A_ori: Adjacency matrix of G_ori, which is the perturbed T, with 0
%            diagonal. A_ori = A_adj_ori - diag(diag(A_adj_ori));



                
