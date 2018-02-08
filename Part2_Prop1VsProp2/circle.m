function circle(x,y,r,numPoint,plotOpt)
% CIRCLE(x,y,r) plots a circle of radius r around point (x,y).
% 
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.01 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)
 
if nargin==3
    plotOpt='.';
    numPoint=51;
elseif nargin==4
    plotOpt='.';
end
 
x=x(:);
y=y(:);
r=r(:);
numCircle=length(x);
theta=linspace(0,2*pi,numPoint)';
xp=reshape(repmat(r,1,numPoint),numCircle*numPoint,1).*repmat(cos(theta),numCircle,1);
yp=reshape(repmat(r,1,numPoint),numCircle*numPoint,1).*repmat(sin(theta),numCircle,1);
plot(repmat(x,numPoint,1)+xp,repmat(y,numPoint,1)+yp,plotOpt);
end