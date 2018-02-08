%% Run 'Main.m' and get 'DataBfFinalPlotRoot39.mat'
%  


close all
clc
clear
numBubbles = 50;
no_region = numBubbles;
T = 60;
delta = 0.2;

stp1 = 3;
[posx,posy,overlapDist] = bubbleSimulatorPublish(numBubbles,stp1);

% load BubbleResult50.mat


[A_adj_ori,A_ori] = GraphTreeGeneration(numBubbles,posx,posy,overlapDist);

stp3=54;
[G_dot_slt,G_slt,n_max,n_slt] = MFDGeneration(no_region,stp3);

stp4=30;
[A_adj,AdjT,Asys,Bsys,L]=Asys_Creation(T,stp4,no_region,A_ori,A_adj_ori,G_dot_slt,G_slt,n_max,n_slt);

plotPlanarTree(no_region,posx,posy,AdjT);



[Kv_r]=DTContTree(T,delta,no_region,A_adj,Asys,Bsys,L);

n0 = n_slt;
root = 39;
sim('LTI_DT_60s_sample')

x_L2normAlongAllPaths(root,AdjT,no_region,x_k_60SampleTime)

