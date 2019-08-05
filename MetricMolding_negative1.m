clear all
a=1; % Edge Length. Keep this as 1. 
%a1=[sqrt(3)*a,0*a];
%a2=[a*sqrt(3)/2, a*3/2];
a1=[a,0];
a2=[0,a];
N=2;
R=5*sqrt(2);
cnt=1
for i= -N:N
    for j= -N:N
%             if  ((a1(1)*i+a2(1)*j)^2 + (a1(2)*i+a2(2)*j)^2 <= (R/1.2)^2 )
                x1(cnt)= a1(1)*i+a2(1)*j;
                y1(cnt)= a1(2)*i+a2(2)*j;
%             end
            
%             if ((a1(1)*i+a2(1)*j)^2+(a1(2)*i+a2(2)*j+a)^2 <= (R/1.2)^2 )
                cnt=cnt+1;
                x1(cnt)= a1(1)*i+a2(1)*j;
                y1(cnt)= a1(2)*i+a2(2)*j+a;
                cnt=cnt+1;
%             end            
    end
end

q=[x1' y1'];
[idx,D]= rangesearch(q,q,a+.1);
%     scatter(x(:),y(:),'.')
%         axis equal
C={};
i=1;

dx=a; % don't change
dy=a; % don't change
scale =dx/2;

syms func1(a,b);
% f(a,b) = sqrt(R^2-(a^2+b^2)); % sphere

f(a,b) = (a^2/R^2-b^2/R^2+1); % negativeCurvature

for cnt_x=1:1:length(x1)
    
    x=x1(cnt_x);
    y=y1(cnt_x);
    
    hx_t=diff(f,a);%-4*exp((-x^2-y^2)/alpha)*x/alpha;
    hy_t=diff(f,b);%-4*exp((-x^2-y^2)/alpha)*y/alpha;
    hx=eval(hx_t(x,y));%x/(sqrt(R^2-(x^2+y^2)));%
    
    hy=eval(hy_t(x,y));%y/(sqrt(R^2-(x^2+y^2)));%
    hxy=hx*hy;
    
    Q=[1+hx^2,hxy/2;...
        hxy/2 1+hy^2];
    [ VectMat, D]=eig(Q);
    
    C{i}={x,y,VectMat,D};
    i=i+1;
end


%% Draw the function
PlotShow=1;
[X,Y] = meshgrid(x1,y1);
Z = eval(f(X,Y));
if PlotShow==1
    surf(X,Y,Z,'FaceAlpha',0.1)
    hold on
end

%%
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
%%
%   k=30
%   cFactor=sqrt(0.5);
    P1= C{k}{3}(:,2)*(E_1(k)/2);
    P2=-1*(C{k}{3}(:,2))*(E_1(k)/2);
    Q1= 1*(C{k}{3}(:,1))*(E_2(k)/2);
    Q2=-1*(C{k}{3}(:,1))*(E_2(k)/2);
    x=C{k}{1}*2*cFactor;
    y=C{k}{2}*2*cFactor;
    coord_X =[ P1(1)+x Q1(1)+x  P2(1)+x  Q2(1)+x ];
    coord_Y= [ P1(2)+y Q1(2)+y  P2(2)+y  Q2(2)+y ];
%   angle = atan2(P2(2),P1(1));
    draw_rectangle([x y],2*cFactor,2*cFactor,0,'w',0.1)
    
    draw_lines([x y],cFactor,E_1(k),E_2(k), eigVect1(:,k), eigVect2(:,k),'k',0.1,0.33) % (Center, cFactor , e1Val , e2Val , E1Vec,E2Vec,rgb,delta,gamma)
    C{i}={x,y,VectMat,D};
    i=i+1;
end
axis equal
axis off
%%
figure
hold on
scale=100;%E_max;
ScalePlot=50
for k=1:length(idx)-1
    k
    P1= C{k}{3}(:,2)*(E_1(k)/2)*scale;
    P2=-1*(C{k}{3}(:,2))*(E_1(k)/2)*scale;
    Q1= 1*(C{k}{3}(:,1))*(E_2(k)/2)*scale;
    Q2=-1*(C{k}{3}(:,1))*(E_2(k)/2)*scale;
    x=C{k}{1}*2*sqrt(c_sq)*ScalePlot;
    y=C{k}{2}*2*sqrt(c_sq)*ScalePlot;
    coord_X =[P1(1)+x Q1(1)+x  P2(1)+x  Q2(1)+x  ];
    coord_Y= [P1(2)+y Q1(2)+y  P2(2)+y  Q2(2)+y ];
%     C{i}={x,y,tmp,D};
    i=i+1;
    plot ([coord_X(1),coord_X(3)], [coord_Y(1),coord_Y(3)],'k');
    plot ([coord_X(2),coord_X(4)], [coord_Y(2),coord_Y(4)],'k');
end

colormap(gray);
%axis([-N/2 N/2 -N/2 N/2])
axis equal
axis off
