function [posx,posy,overlapDist,N] = bubbleSimulatorPublish(numBubbles,stp1)
% parameters
% numBubbles=10;     % number of bubbles
minRadius=200;       % minimum radius
maxRadius=300;     % maximum radius
k=50;               % spring constant
c=2;                % damping constant
cg=1;               % damping w.r.t. ground
g=50;               % gravity
Ts=0.003;           % simulation sampling time
rand('seed',stp1)
 
bubbles.radius=rand(numBubbles,1)*(maxRadius-minRadius)+minRadius;
bubbles.mass=0.1*ones(numBubbles,1);
bubbles.pos=(rand(numBubbles,2)-0.5)*sqrt(2)*sqrt(numBubbles)*maxRadius;
bubbles.vel=zeros(numBubbles,2);
 
% offline computation
rr=repmat(bubbles.radius,1,numBubbles);
bubbles.sumRadius=(rr+rr');
 
t=0;
% circle(bubbles.pos(:,1),bubbles.pos(:,2),bubbles.radius);
% axis equal;drawnow;
% simulation starts
tic
while t<10 
    t=t+Ts;
    [posx,posy,overlapDist,N,bubbles]=updatePosition(bubbles,numBubbles,Ts,k,c,cg,g,t);
%     if ~mod(round(t/Ts),10)
%         circle(bubbles.pos(:,1),bubbles.pos(:,2),bubbles.radius);
%         axis equal;drawnow;
%     end
end
xx = toc

end

function [posx,posy,overlapDist,N,bubbles]=updatePosition(bubbles,numBubbles,Ts,k,c,cg,g,t)
% consider elasticity, damping relative to ground and attraction to the centre
posx=bubbles.pos(:,1);
posy=bubbles.pos(:,2);
relPosx=repmat(posx',numBubbles,1)-repmat(posx,1,numBubbles);
relPosy=repmat(posy',numBubbles,1)-repmat(posy,1,numBubbles);
relDist=sqrt(relPosx.^2+relPosy.^2);
overlapDist=bubbles.sumRadius-relDist; 
overlapDist(1:numBubbles+1:numBubbles*numBubbles)=-inf;
velx=bubbles.vel(:,1);
vely=bubbles.vel(:,2);
relVelx=repmat(velx',numBubbles,1)-repmat(velx,1,numBubbles);
relVely=repmat(vely',numBubbles,1)-repmat(vely,1,numBubbles);
relSpeed=sqrt(relVelx.^2+relVely.^2);
 
% calculate forces and acceleration
N=overlapDist>0; 
Ns=N & (relSpeed>0);
F=zeros(numBubbles,2);
C=zeros(numBubbles,2);
for i=1:numBubbles
    Ni=N(i,:);
    overlapDistNi=overlapDist(i,Ni);
    relDistNi=relDist(i,Ni);
    normedF=(k*overlapDistNi)./relDistNi;
    F(i,:)=normedF*[-relPosx(i,Ni)' -relPosy(i,Ni)'];
     
    Nsi=Ns(i);
    relSpeedNsi=relSpeed(i,Nsi);
    normedC=(c*relSpeedNsi)./relSpeedNsi;
    C(i,:)=normedC*[relVelx(i,Nsi)' relVely(i,Nsi)'];
end
Cg=-cg*bubbles.vel;
G=-g*bubbles.pos*0.05;
accel=(F+C+Cg)./repmat(bubbles.mass,1,2)+G;
 
% update states
bubbles.pos=Ts*bubbles.vel+bubbles.pos;
bubbles.vel=Ts*accel+bubbles.vel;

end