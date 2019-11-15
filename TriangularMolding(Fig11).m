%% Introduction
%{
Step 1: Input a target function and plot it.
Step 2: Draw the triangular Grid on which we draw the precursor
pattern'
Step 3: Draw the Precursor pattern.
We prepare a grid and draw a target 3D function on our grid that we 
wish to generate.
Then we create the 2-D pattern that will give us the target 3D pattern on 
contracting out the black regions (painted yellow here).
%}

%% Step 1 : Input the target Function and plot on a 3d Grid
clear all

num=5;
iVal=num;
jVal=num;
a=1;
x0=0;
y0=0;
CoordMat=zeros(iVal,jVal,2);
for i=1:1:iVal
    for j=1:1:jVal
        if rem(j,2)==0
            CoordMat(i,j,1)=x0+i*a;
        elseif rem(j,2)==1
            CoordMat(i,j,1)=x0+i*a+a/2;
%         elseif rem(j,3)==2
%             CoordMat(i,j,1)=x0+i*a+a;
        end
        CoordMat(i,j,2)=y0+sqrt(3)/2*a*j;
    end
end
% figure();
% scatter(reshape(CoordMat(:,:,1),[num^2,1]),reshape(CoordMat(:,:,2),[num^2,1]),'.');
a_=CoordMat(uint8(num/2),uint8(num/2),1); b_=CoordMat(uint8(num/2),uint8(num/2),2); R_=3; %Parameters of circle
 
X_vertices=reshape(CoordMat(:,:,1),[num^2,1]);
Y_vertices=reshape(CoordMat(:,:,2),[num^2,1]);

syms funcz(a,b);
R=3;                      
% f(a,b) = (a^2/R^2-b^2/R^2+1);   % saddle or negativeCurvature
funcz(a,b) = sqrt(R^2-((a-3).^2+(b-3).^2));     % sphere or Positive Curvature

Z_vertices=real(eval(funcz(X_vertices,Y_vertices)));
figure();
scatter3(X_vertices,Y_vertices,Z_vertices,'filled')
title('Target Function in 3d Scatter on the grid')
%% Step 2:Draw Triangular Lattice
figure
X_center=[];Y_center=[];X_vert=[];Y_vert=[];
m=1; %n=1;
for j = 1:1:jVal-1
    for i=1:1:iVal-1
        if rem(j,2)==1
            % Triangle 1 is formed by (i,j) and (i+1,j) and (i,j+1)
            % Find final size of the triangle
            A1=CoordMat(i,j,1); A2=CoordMat(i,j,2);
            B1=CoordMat(i,j+1,1); B2=CoordMat(i,j+1,2);
            C1=CoordMat(i+1,j+1,1); C2=CoordMat(i+1,j+1,2);
            plot([A1,B1,C1,A1],[A2,B2,C2,A2],'r')
            X_center(m)=(A1+B1+C1)/3; Y_center(m)=(A2+B2+C2)/3;
            X_vert(m,:)=[A1;B1;C1];
            Y_vert(m,:)=[A2;B2;C2];
            m=m+1;
            
            A1=CoordMat(i,j,1); A2=CoordMat(i,j,2);
            B1=CoordMat(i+1,j,1); B2=CoordMat(i+1,j,2);
            C1=CoordMat(i+1,j+1,1); C2=CoordMat(i+1,j+1,2);
            plot([A1,B1,C1,A1],[A2,B2,C2,A2],'r')
            X_center(m)=(A1+B1+C1)/3; Y_center(m)=(A2+B2+C2)/3;
            X_vert(m,:)=[A1;B1;C1];
            Y_vert(m,:)=[A2;B2;C2];
            m=m+1;
 
        else
            % Triangle 1 is formed by (i,j) and (i+1,j) and (i,j+1)
            % Find final size of the triangle
            A1=CoordMat(i,j,1); A2=CoordMat(i,j,2);
            B1=CoordMat(i,j+1,1); B2=CoordMat(i,j+1,2);
            C1=CoordMat(i+1,j,1); C2=CoordMat(i+1,j,2);
            plot([A1,B1,C1,A1],[A2,B2,C2,A2],'r')
            X_center(m)=(A1+B1+C1)/3; Y_center(m)=(A2+B2+C2)/3;
            X_vert(m,:)=[A1;B1;C1];
            Y_vert(m,:)=[A2;B2;C2];
            m=m+1;
            A1=CoordMat(i+1,j+1,1); A2=CoordMat(i+1,j+1,2);
            B1=CoordMat(i,j+1,1); B2=CoordMat(i,j+1,2);
            C1=CoordMat(i+1,j,1); C2=CoordMat(i+1,j,2);
            plot([A1,B1,C1,A1],[A2,B2,C2,A2],'r')
            X_center(m)=(A1+B1+C1)/3; Y_center(m)=(A2+B2+C2)/3;
            X_vert(m,:)=[A1;B1;C1];
            Y_vert(m,:)=[A2;B2;C2];
            m=m+1;
        end
    end
end
X_vert=X_vert';
Y_vert=Y_vert';
 
Z=arrayfun(@(x1,x2) funcz(x1,x2),X_vert,Y_vert);
% scatter3(X_vert(:),Y_vert(:),Z(:))
% figure
hold on
 
axis equal
scatter(X_center(:),Y_center(:),'r.')
hexagons = plot([X_vert;X_vert(1,:)],[Y_vert;Y_vert(1,:)],'r-');
axis off
 
title('3D Grid to draw the pattern on')

%% Step 3:Generate the precursor pattern
figure
Xc=[X_center;X_center;X_center]; % Center Coordinates
Yc=[Y_center;Y_center;Y_center];
CenterLengths=zeros(size(Xc));
for p=1:1:numel(Xc)
    CenterLengths(p)=norm([X_vert(p)-Xc(p),Y_vert(p)-Yc(p)]);
end
NewCenterLengths=zeros(size(Xc));
 
Z_vert=ones([3,length(X_center)]);
 
Zc=[funcz(X_center,Y_center);funcz(X_center,Y_center);funcz(X_center,Y_center)];
for p=1:1:numel(Xc)
    Z_vert(p)=funcz(X_vert(p),Y_vert(p));
    NewCenterLengths(p)=norm([X_vert(p)-Xc(p),Y_vert(p)-Yc(p),Z_vert(p)-Zc(p)]);
end
 
xEdge=[X_vert(1,:)-X_vert(2,:);X_vert(2,:)-X_vert(3,:);X_vert(3,:)-X_vert(1,:)];
yEdge=[Y_vert(1,:)-Y_vert(2,:);Y_vert(2,:)-Y_vert(3,:);Y_vert(3,:)-Y_vert(1,:)];
EdgeLengths=ones(3,length(X_center));
for p=1:1:numel(Xc)
    EdgeLengths(p)=norm([xEdge(p),yEdge(p)]);
end
 
zEdge=[Z_vert(1,:)-Z_vert(2,:);Z_vert(2,:)-Z_vert(3,:);Z_vert(3,:)-Z_vert(1,:)];
NewEdgeLengths=ones(3,length(X_center));
for p=1:1:numel(Xc)
    NewEdgeLengths(p)=norm([xEdge(p),yEdge(p),zEdge(p)]);
end
Cmat=[NewEdgeLengths./EdgeLengths;NewCenterLengths./CenterLengths];
CmatSort=sort(reshape(Cmat,[],1));
CFactor=CmatSort(end);  % Selecting the maximum value.
% We can sometimes set it to 
% CmatSort(end-1) % as well, so as to reduce the black region.
 

 
%plot
% figure
% scatter(X_vert(:,i),Y_vert(:,i))
hold on
% We create the new matrix with the new side length
% figure
for i=1:1:length(Xc)
    Cijk=[X_center(i),Y_center(i)];
    Ai=[X_vert(1,i),Y_vert(1,i)];
    Aj=[X_vert(2,i),Y_vert(2,i)];
    Ak=[X_vert(3,i),Y_vert(3,i)];
    
    Dij=(Ai+Aj)/2-(1-Cmat(1,i)/CFactor)/2*(Aj-Ai);
    Dji=(Ai+Aj)/2+(1-Cmat(1,i)/CFactor)/2*(Aj-Ai);
 
    Djk=(Aj+Ak)/2-(1-Cmat(2,i)/CFactor)/2*(Ak-Aj);
    Dkj=(Aj+Ak)/2+(1-Cmat(2,i)/CFactor)/2*(Ak-Aj);
 
    Dik=(Ai+Ak)/2-(1-Cmat(3,i)/CFactor)/2*(Ak-Ai);
    Dki=(Ai+Ak)/2+(1-Cmat(3,i)/CFactor)/2*(Ak-Ai);
    
    
    Ei=Cijk+(1-Cmat(4,i)/CFactor)*(Ai-Cijk);
    Ej=Cijk+(1-Cmat(5,i)/CFactor)*(Aj-Cijk);
    Ek=Cijk+(1-Cmat(6,i)/CFactor)*(Ak-Cijk);
 
    X=[Djk(1),Dkj(1),Ek(1),Dki(1),Dik(1),Ei(1),Dij(1),Dji(1),Ej(1),Djk(1)];
    Y=[Djk(2),Dkj(2),Ek(2),Dki(2),Dik(2),Ei(2),Dij(2),Dji(2),Ej(2),Djk(2)];
    fill(X,Y,'y')
end
hold off
axis equal
axis off
% We find the final length of the three sides and 3 radii
% We select the points of interest and draw the line.
%Fill in the colors? Optional
'Yellow portion is to be painted black. Rest to be painted white.'
title('Required precursor pattern')
