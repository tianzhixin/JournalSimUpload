function [dtaA,Dtax] = CalxDiff(Aori,Atree,x_tree,x_scp_ori_fst)

%% Calculate ||deltaA||_{L,1}
dtaA = sum(sum(abs(Aori - Atree)));

%% Get Scope Data
xTree = x_tree.signals.values;
xOrik = x_scp_ori_fst.signals.values;
Deltaxk = xOrik - xTree;
Dtax = sqrt(sum(sum(Deltaxk.^2)));

end
