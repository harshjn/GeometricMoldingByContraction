% Surface Revolution of a function will give us the target 3d shape that we wish to generate via tailoring
% We can choose to split the patterns towards the end in order to reduce
% the effect of wrinkling at ends with the boolean SPLIT.

clear all; clc; close all
deltaX=0.01;  % Set this to smaller value for a more continous pattern.
xRangeStart=deltaX; xRangeEnd=1.4;
x_Mat=xRangeStart:deltaX:xRangeEnd;
A=1;%B=1;
% Define the function which is to be revolved around y-axis.
func=@(x) A*(x).^2;% Function to be revolved to obtain the surface

% Use func2 as well in case of piecewise fuctions
% syms func2(x);
% func2(x) = (piecewise(0<x <0.3, func(0.3), 0.3 <= x <= xRangeEnd, func(x), 1));

y_Mat=func(x_Mat);
plot(x_Mat,y_Mat,'.')
title('surface to be revolved')

SPLIT=0;

%% Plot the fraction of circumference to be removed as a function of radius.
figure;
i=length(x_Mat);
Per_sub=zeros(i,1);
Frac=zeros(i,1);
r_Mat=zeros(i,1);

for j=1:1:i
    dx=0.01;  % Can be set to lower value such as 0.001 for more accurate results
    Xmat=0:dx:x_Mat(j);
%     Ymat=-A*Xmat.^2+B*Xmat.^4;
    Ymat=(func(Xmat));
    r=arclength(Xmat,Ymat,'sp');
    r_Mat(j)=r;
    Per_sub(j)=2*pi*(r-x_Mat(j));
    Frac(j)=2*pi*(r-x_Mat(j))./(2*pi*r);
end
plot(r_Mat,Frac)
title('Fraction of perimeter to be removed at radius r')
%% Distribute in n-Petals with optional SPLIT
% We divide the fraction to be removed into petals.
figure
% Contraction factor gamma given as InitalLength/ContractedLength
gamma=0.5;         
n_=8; % Number of petals
alpha=1/gamma;  % Could be replaced by a matrix instead of a constant value
alphaMat=alpha*ones(length(r_Mat),1);


cycleCount=length(r_Mat)

if SPLIT==1
%   Split the alpha for different petals.
    sin_Mat=alpha*sin([0:pi/length(r_Mat):pi-pi/length(r_Mat)]);
    sin_Mat(1:length(r_Mat)/2)=alpha;
    alpha_=sin_Mat'.*ones(length(r_Mat),1);
    Left_Mat=alphaMat-alpha_;
    rMat=r_Mat;
    FracC=alpha_.*Frac;
    FracLeft=Left_Mat.*Frac;
    
    for i=2:1:length(rMat)
        i
    %     %Increase the number of petals in a range for rMat optionally.
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
    
else
    
    FracC=alpha.*Frac;
    rMat=r_Mat;
    for i=2:1:length(rMat)
    i
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

    end     
    end
    grid off
    set(gca,'RTickLabel',[]);
end
