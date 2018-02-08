close all
clc
clear
numBubbles = 50;
no_region = numBubbles;
T = 60;
delta = 0.2;

stp1 = 3;
[posx,posy,overlapDist] = bubbleSimulatorPublish(numBubbles,stp1);

[A_adj_ori,A_ori] = GraphTreeGeneration(numBubbles,posx,posy,overlapDist);

stp3=54;
[G_dot_slt,G_slt,n_max,n_slt] = MFDGeneration(no_region,stp3);

stp4=30;
[Asys_ori,Bsys_ori,L]=Asys_ori_Creation(T,stp4,no_region,A_ori,A_adj_ori,G_dot_slt,G_slt,n_max,n_slt);

NumTree = 15;
deltaA_L1 = zeros(NumTree,1);
DeltaxEgy = deltaA_L1;
disSeedArr = 1:no_region; 

for i = 1:NumTree
    waitbar(i / NumTree)

    if i==1
        [A_adj,AdjT,Asys,Bsys] = Asys_generation(Asys_ori,no_region,G_slt);
    else
        [A_adj,AdjT,Asys,Bsys] = GenerSpanningTsimtoAsys(i,no_region,Asys_ori,G_slt);
    end
    [Kv_r]=DTContTree(T,delta,no_region,A_adj,Asys,Bsys,L);
    n0 = n_slt;
    Atree = eye(no_region)+T*L*(Asys + Bsys*Kv_r)*inv(L);
    Aori = eye(no_region)+T*L*(Asys_ori + Bsys*Kv_r)*inv(L);
    sim('LTI_DT_60s_sample')
    [deltaA_L1(i),DeltaxEgy(i)] = CalxDiff(Aori,Atree,x_tree,x_scp_ori_fst);
end



% plotPlanarTree(no_region,posx,posy,AdjT);







