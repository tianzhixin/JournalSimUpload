%% This code is to generate the path from a specified root to all leaves.
%
%% The flow goes like this:
% Get the level of each node given a root i.
% Generate the adjacency matrix for a directed tree with root i.
% Find all leaves
% Find all paths from root i and save in cell Path

%% Get the level of each node given a root i.
% load AsysBsys50Regions_Tree170824.mat
function Path = PathFinding(AdjT,root,no_region)

[dist] = graphallshortestpaths(AdjT);
Height = max(dist(root,:));
Level = dist(root,:);

%% Generate the adjacency matrix for a directed tree with root i.
AdjTrt = AdjT;
for i = 1:no_region
    for j = 1:no_region
        if AdjTrt(i,j)~=0
            if Level(i) > Level(j)
                AdjTrt(i,j) = 0;
            end
        end
    end
end
% g1 = view(biograph(AdjTrt))

%% Find all leaves

AdjTfull = full(AdjT);
numLeaves = 0;
for i = 1:no_region
    if length(find(AdjTfull(i,:))) == 1
        if i ~= root
            Leaves(numLeaves+1) = i;
            numLeaves = numLeaves+1;
        end
    end
end

% if length(find(AdjTfull(root,:))) == 1
%     numPaths = numLeaves - 1;
% else
%     numPaths = numLeaves;
% end

%% Find all paths  from root i and save in cell Path

Path = cell(numLeaves,1);

for i = 1:numLeaves
    EndNode = Leaves(i);
    fromRoot = graphtraverse(AdjTrt,root);
    toEndNode = graphtraverse(AdjTrt',EndNode);
    h = intersect(fromRoot,toEndNode);
    Path2i = zeros(size(h));
    for j = 1:length(h)
        Path2i(Level(h(j))+1) = h(j);
    end
    Path{i} = Path2i;
end
        
        
    
 %% Test code to show how to find all paths   

% EndNode = 13;
% fromRoot = graphtraverse(AdjTrt,root);
% toEndNode = graphtraverse(AdjTrt',EndNode);
% h = intersect(fromRoot,toEndNode);
% DG2 = AdjTrt(h,h);
% g2 = view(biograph(DG2,cellstr(num2str(h'))))           
