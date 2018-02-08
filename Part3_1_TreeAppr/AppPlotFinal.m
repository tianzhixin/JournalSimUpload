close all

load appTest;
deltaA_L1_all = deltaA_L1;
DeltaxEgy_all = DeltaxEgy;

load appTest6;
deltaA_L1_all = [deltaA_L1_all;deltaA_L1];
DeltaxEgy_all = [DeltaxEgy_all;DeltaxEgy];

figure
% B = regress(DeltaxEgy,deltaA_L1);
% z = B(1)*deltaA_L1;
plot(deltaA_L1_all(2:length(deltaA_L1_all)),DeltaxEgy_all(2:length(DeltaxEgy_all)),'bo','markersize',10,'linewidth',2)
hold on
plot(deltaA_L1_all(1),DeltaxEgy_all(1),'ro','markersize',10,'linewidth',2)
xlabel('$\sum_{i=1}^{m} \sum_{j=1}^{m} |\Delta A_{i,j}|$','Interpreter','LaTex')
ylabel('$\delta$','Interpreter','LaTex')