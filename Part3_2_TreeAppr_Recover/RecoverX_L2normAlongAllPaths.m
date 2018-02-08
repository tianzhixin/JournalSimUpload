
function RecoverX_L2normAlongAllPaths(root,AdjT,no_region,x_scp_ori_fst)

%% delta_nt_nor & delta_nk_nor calculation & plot
% 
% load AsysBsys50Regions_Tree170824.mat
% root = 4;
Path = PathFinding(AdjT,root,no_region);
[dist] = graphallshortestpaths(AdjT);
Level = dist(root,:);


x_60 = x_scp_ori_fst.signals.values;


Egy = zeros(1,no_region);
EgNoLog =  zeros(1,no_region);
for i = 1:no_region
    Egy(i) = log10(sqrt(sum((x_60(:,i).^2))));
    EgyNoLog(i) = sqrt(sum((x_60(:,i).^2)));

end

figure
for i = 1:length(Path)
    pt = Path{i,1};
    plot(Level(pt),Egy(pt),'*-','LineWidth',1.5);
    hold on;
end
xlabel('Level of node in $\mathcal{T}_{39}$','Interpreter','LaTex')
ylabel('$\gamma$','Interpreter','LaTex')


% legend({'Path 1','Path 2','Path 3','Path 4','Path 5',...
%     'Path 6','Path 7','Path 8','Path 9','Path 10'},'Location','NorthEast');
% 
% legend({'Path 11','Path 12','Path 13','Path 14','Path 15',...
%     'Path 16','Path 17','Path 18','Path 19','Path 20'},'Location','NorthEast');


% legend({'Path 1','Path 2','Path 3','Path 4','Path 5',...
%     'Path 6','Path 7','Path 10','Path 9','Path 10','Path 11','Path 12','Path 13','Path 14','Path 15',...
%     'Path 16','Path 17','Path 18','Path 19','Path 20'},'Location','EastOutside')    

% figure
% for i = 1:length(Path)
%     pt = Path{i,1};
%     plot(Level(pt),EgyNoLog(pt),'*-','LineWidth',1.5);
%     hold on;
% end
end
%% --------Draw delta_n(k) in tree unstable case--------------- 
