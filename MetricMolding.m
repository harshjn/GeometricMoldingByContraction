
%% Introduction
%{
Uses the functions draw_rectangle and draw_lines and RescalePoly.

We prepare a grid and draw a target 3D function on our grid that we 
wish to generate.
Then we create the 2-D pattern that will give us the target 3D pattern on contracting out the black regions.

The thickness of the black boundaries can be controlled by factor delta.
The Gamma factor for contraction property of the material can also be controlled. 
gamma = (final length after contraction)/(initial length before contraction)

%}

%% Step 1
% Prepare a grid, plot our target function on it.
clear all
a=1; % Edge Length. Keep this as 1. Image can be scaled later.
a1=[a,0]; % basis vectors for our grid
a2=[0,a];

TargetSize=5; % If we need a 5x5 grid
N=floor(TargetSize/2);

count=1

for i= -N:N
    for j= -N:N
%             if  ((a1(1)*i+a2(1)*j)^2 + (a1(2)*i+a2(2)*j)^2 <= (R/1.2)^2 )
                x1(count)= a1(1)*i+a2(1)*j;
                y1(count)= a1(2)*i+a2(2)*j;
%             end
            
%             if ((a1(1)*i+a2(1)*j)^2+(a1(2)*i+a2(2)*j+a)^2 <= (R/1.2)^2 )
                count=count+1;
                x1(count)= a1(1)*i+a2(1)*j;
                y1(count)= a1(2)*i+a2(2)*j+a;
                count=count+1;
%             end            
    end
end

q=[x1' y1'];
[idx,D]= rangesearch(q,q,a+.1);
%     scatter(x(:),y(:),'.')
%         axis equal
C={};
i=1;

dx=a; dy=a;
scale =dx/2;

syms func1(a,b);
R=5*sqrt(2);
f(a,b) = (a^2/R^2-b^2/R^2+1); % negativeCurvature
% f(a,b) = sqrt(R^2-(a^2+b^2)); % sphere

PlotShow=1;
[X,Y] = meshgrid(x1,y1);
Z = eval(f(X,Y));
if PlotShow==1
    surf(X,Y,Z,'FaceAlpha',0.1)
    hold on
end
title('Target 3D function on the 2D grid')
% scatter(x1,y1) % Plot the grid
view(-121,28)
%% Generating the 2D pattern that gives target 3D pattern on heating.
for count_x=1:1:length(x1)
    x=x1(count_x);       y=y1(count_x);
    hx_t=diff(f,a);      hy_t=diff(f,b);
    hx=eval(hx_t(x,y));  hy=eval(hy_t(x,y));
    hxy=hx*hy;
    Q=[1+hx^2,hxy/2;...
        hxy/2 1+hy^2];
    [ VectMat, D]=eig(Q);
    C{i}={x,y,VectMat,D};
    i=i+1;
end

figure
hold on
eigVect1=zeros(2,length(idx));
eigVect2=zeros(2,length(idx));
for k= 1: length(idx)
    eigen_1(k)= C{k}{4}(1);
    eigen_2(k)= C{k}{4}(4);
    eigVect1(:,k)=C{k}{3}(:,1);
    eigVect2(:,k)=C{k}{3}(:,2);
    
end
E_max=max(max(eigen_1),max(eigen_2)); % this value is equal to c
c_sq=1.0*E_max;
cFactor=sqrt(c_sq);
E_1=eigen_1;
E_2=eigen_2;
for k=1:length(idx)-1
    P1= C{k}{3}(:,2)*(E_1(k)/2);
    P2=-1*(C{k}{3}(:,2))*(E_1(k)/2);
    Q1= 1*(C{k}{3}(:,1))*(E_2(k)/2);
    Q2=-1*(C{k}{3}(:,1))*(E_2(k)/2);
    x=C{k}{1}*2*cFactor;
    y=C{k}{2}*2*cFactor;
    coord_X =[ P1(1)+x Q1(1)+x  P2(1)+x  Q2(1)+x ];
    coord_Y= [ P1(2)+y Q1(2)+y  P2(2)+y  Q2(2)+y ];
    draw_rectangle([x y],2*cFactor,2*cFactor,0,'w',0.1)
    % We input contraction factor gamma as a fraction from 0 to 1.                                              
    % gamma=Final/Initial Length
    % We also input the fraction of black boundary surrounding the regions.
    draw_lines([x y],cFactor,E_1(k),E_2(k), eigVect1(:,k), eigVect2(:,k),'k',0.1,0.33) 
            % (Center, cFactor , e1Val , e2Val , E1Vec,E2Vec,rgb,delta,gamma)
    C{i}={x,y,VectMat,D};
    i=i+1;
end
axis equal
axis off
