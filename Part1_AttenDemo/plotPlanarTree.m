function plotPlanarTree(no_region,posx,posy,AdjT)

figure
for i=1:no_region
    for j=1:i
        if AdjT(i,j)                                           %N=overlapDist>0; 
        hold on, plot([posx(i),posx(j)],[posy(i),posy(j)],'LineWidth',1.5,'color',[0 0.5 1])          % If circle i and circle j overlap, then draw an edge connecting them.
        end
    end
end

plot(posx,posy,'.','Markersize',30,'color',[0 0 1])

for i=1:no_region
    txt{i}=num2str(i); 
end



for i = 1:no_region
%     if i==47
%         text(posx(i)-180,posy(i)-80,txt{i},'FontSize',13)
%     else
%         text(posx(i)+40,posy(i)+40,txt{i},'FontSize',13)
%     end
        text(posx(i)+40,posy(i)+40,txt{i},'FontSize',13)
end

end