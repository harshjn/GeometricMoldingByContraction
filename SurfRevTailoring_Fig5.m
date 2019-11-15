% Surface Revolution of a function will give us the target 3d shape that we wish to generate via tailoring
clear all; clc; close all
deltaX=0.001;
xRangeStart=deltaX; xRangeEnd=1.4;
x_Mat=xRangeStart:deltaX:xRangeEnd;
A=1;%B=1;

% Define the function which is to be revolved around y-axis.
func=@(x) -A*(x-0.85).^3;% Function to be revolved to obtain the surface

% Use func2 as well in case of piecewise fuctions
% syms func2(x);
% func2(x) = (piecewise(0<x <0.3, func(0.3), 0.3 <= x <= xRangeEnd, func(x), 1));

y_Mat=func(x_Mat);
plot(x_Mat,y_Mat,'.')
title('surface to be revolved')

%% Plot the fraction of circumference to be removed as a function of radius.
i=length(x_Mat);
Per_sub=zeros(i,1);
Frac=zeros(i,1);
r_Mat=zeros(i,1);

for j=1:1:i
    dx=0.001;
    Xmat=0:dx:x_Mat(j);
%     Ymat=-A*Xmat.^2+B*Xmat.^4;
    Ymat=eval(func(Xmat));
    r=arclength(Xmat,Ymat,'sp');
    r_Mat(j)=r;
    Per_sub(j)=2*pi*(r-x_Mat(j));
    Frac(j)=2*pi*(r-x_Mat(j))./(2*pi*r);
end
plot(r_Mat,Frac)

%% Distribute in n-Petals
figure
alpha=2; % Contraction factor given as InitalLength/ContractedLength
          % Could be replaced by a matrix instead of a constant value
alphaMat=alpha*ones(length(r_Mat),1);
sin_Mat=alpha*sin([0:pi/length(r_Mat):pi-pi/length(r_Mat)]);
% sin_Mat(1:length(r_Mat)/2)=alpha;
alpha_=sin_Mat'.*ones(length(r_Mat),1);
Left_Mat=alphaMat-alpha_;
rMat=r_Mat;
FracC=alpha_.*Frac;
n_=8;
FracLeft=Left_Mat.*Frac;
for i=2:1:length(rMat)
%     %Generate the Polar angle vector containing information about sector location and angle
%     if i>length(rMat)/2
%         n=2*n_;
%     else
%         n=n_;
%     end
    n=n_;
    ThetaMat=0:2*pi/n:2*pi-2*pi/n;
    for j=1:1:n
        Theta = ThetaMat(j)-FracC(i)*2*pi/(2*n):0.0001:ThetaMat(j)+FracC(i)*2*pi/(2*n);  
        %Generate the Radius vector
        R = rMat(i);
        %Create a grid from angle and Radius
        [ThetaG,RG] = meshgrid(Theta,R);
        %Create X,Y matrices calculated on grid.
        X = RG.*cos(ThetaG);
        Y = RG.*sin(ThetaG);
        %Calculate the function
        Z = RG;
        %Surf plot
        polarplot(ThetaG,RG,'.k','MarkerSize',3);
        hold on
        
        if(i>length(rMat)/2)
            Theta = ThetaMat(j)+(2*pi/(2*n))-FracLeft(i)*2*pi/(2*n):0.0001:...
                ThetaMat(j)+(2*pi/(2*n))+FracLeft(i)*2*pi/(2*n);  
            %Generate the Radius vector
            R = rMat(i);
            %Create a grid from angle and Radius
            [ThetaG,RG] = meshgrid(Theta,R);
            %Create X,Y matrices calculated on grid.
            X = RG.*cos(ThetaG);
            Y = RG.*sin(ThetaG);
            %Calculate the function
            Z = RG;
            %Surf plot
            polarplot(ThetaG,RG,'.k','MarkerSize',3);
            hold on
        end
    end     
end
grid off
set(gca,'RTickLabel',[]);

