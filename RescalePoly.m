function [A,B]= RescalePoly(InCoord,CentCoord,Factor)

OutCoord=zeros(size(InCoord));
for k=1:1: size(InCoord,1)
OutCoord(k,:)=CentCoord+(InCoord(k,:)-CentCoord)*Factor;
end

A=OutCoord(:,1);B=OutCoord(:,2);
end
