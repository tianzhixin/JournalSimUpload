close all
clc
clear

NwSize = [20 30 40 50 10];  % Network size
NumT = 10;                  % For each network size, generate NumT different tree structures.

T = 60;                     % Sampling time
delta = 0.2;                % sum(Q_i,j) < 1- delta


tProp1 = zeros(NumT,length(NwSize));        % Matrix recording time consumption for P1 with row index of structure and column index of network size
tProp2 = tProp1;                            % Matrix recording time consumption for P2

RunTimes1 = 10;                              % Run 'RunTimes' times and get the mean value
RunTimes2 = 10;

t1 = zeros(length(NwSize),NumT,RunTimes1);
t2 = zeros(length(NwSize),NumT,RunTimes2);
stp1=82;


for sizeIdx = 1:length(NwSize)
    numBubbles = NwSize(sizeIdx);
    no_region = numBubbles;
    
    for stp3 = 1:NumT
        close all;
        waitbar((NumT*(sizeIdx-1)+stp3)/(NumT*length(NwSize)))

        [posx,posy,overlapDist] = bubbleSimulatorPublish(numBubbles,stp1);
        
        [A_adj_ori,A_ori] = GraphTreeGeneration(numBubbles,posx,posy,overlapDist);
        
        [G_dot_slt,G_slt,n_max,n_slt] = MFDGeneration(no_region,stp3);
        
        stp4=1;
        [A_adj,AdjT,Asys,Bsys,L]=Asys_Creation(T,stp4,no_region,A_ori,A_adj_ori,G_dot_slt,G_slt,n_max,n_slt);

        [Kv_r]=DTContTree(T,delta,no_region,A_adj,Asys,Bsys,L);
        
        %% ----------Time consumption of Prop2----------------

        Qv = L*(eye(no_region)+T*(Asys+Bsys*Kv_r))*inv(L);

        

        %% ----------Time consumption of Prop1----------------

        h = @() Prop1time(Qv,no_region,T,A_adj);
    
        
        suc = Prop1time(Qv,no_region,T,A_adj);
        if suc == 1
            for i = 1:RunTimes1
                t1(sizeIdx,stp3,i) = timeit(h);
            end
            tProp1(stp3,sizeIdx) = mean(t1(sizeIdx,stp3,:));
        end
        
        
        if suc == 1
            g = @() Prop2time(Qv,no_region);
            for i = 1:RunTimes2
                t2(sizeIdx,stp3,i) = timeit(g);
            end
            tProp2(stp3,sizeIdx) = mean(t2(sizeIdx,stp3,:));
        end
        
    end
end

[r,c] = size(tProp1);
tProp1m = zeros(1,c);
for i = 1:c
    tProp1m(i) = sum(tProp1(:,i))/r;
end

tProp2m = zeros(1,c);
for i = 1:c
    tProp2m(i) = sum(tProp2(:,i))/r;
end

% t1 = zeros(length(NwSize),NumT,RunTimes1);

for i = 1:length(NwSize)
    t1Mapdev = t1(i,:,:) - tProp1m(i)*ones(1,NumT,RunTimes1);
    t1dev(i) = sqrt(sum(sum(t1Mapdev.^2))/(NumT*RunTimes1-1));
    t2Mapdev = t2(i,:,:) - tProp2m(i)*ones(1,NumT,RunTimes2);
    t2dev(i) = sqrt(sum(sum(t2Mapdev.^2))/(NumT*RunTimes2-1));
end

% load BubbleResult50.mat






% plotPlanarTree(no_region,posx,posy,AdjT);






