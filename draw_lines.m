function[]= draw_lines(Center, cFactor , e1Val , e2Val , E1Vec,E2Vec,rgb,delta,gamma)

Factor=(1-delta/2);
% E1Vec=[1 2]; E2Vec=[-4 3];
Edge=cFactor;
E1Vec=E1Vec/norm(E1Vec); E2Vec=E2Vec/norm(E2Vec);
% e1Val=1;e2Val=1;
PE1=[-E1Vec(2) E1Vec(1)];
PE1=PE1/norm(PE1); % Perpendicular to E1
PE2=[-E2Vec(2) E2Vec(1)];
PE2=PE2/norm(PE2);
angle1=atan2(E1Vec(2),E1Vec(1)); angle2=atan2(E2Vec(2),E2Vec(1));
alpha= Edge*(1-sqrt(e1Val)/cFactor)/(1-gamma); % gamma=final/initial length
beta = Edge*(1-sqrt(e2Val)/cFactor)/(1-gamma);
%Edge=Delta-delta;
% Edge=3;
% for k = 1%0:dx:GridEdge
C11=Center+alpha/2*PE1;
C12=Center-alpha/2*PE1;
C21=Center+beta/2*PE2;
C22=Center-beta/2*PE2;
D1=C11+E1Vec'*Edge*(2);D2=C12+E1Vec'*Edge*2;
B1=C11-E1Vec'*Edge*(2);B2=C12-E1Vec'*Edge*2;
%     fill([B1(1),B2(1),D2(1),D1(1)],[B1(2),B2(2),D2(2),D1(2)],'k')
%     hold on
A1=C21-E2Vec'*Edge*(2);A2=C22-E2Vec'*Edge*2;
C1=C21+E2Vec'*Edge*(2);C2=C22+E2Vec'*Edge*2;
%     fill([A1(1),A2(1),C2(1),C1(1)],[A1(2),A2(2),C2(2),C1(2)],'k')
%     rectangle('Position',[Center(1)-Edge Center(2)-Edge 2*Edge 2*Edge],'EdgeColor','r')
%     hold off
%%
L=Edge;H=Edge; %Center=Center;
xlimit = [Center(1)-L Center(1)+L];
ylimit = [Center(2)-H  Center(2)+H];
xbox = xlimit([1 1 2 2 1]);
ybox = ylimit([1 2 2 1 1]);

%     mapshow(xbox,ybox,'DisplayType','polygon','LineStyle','none')
hold on
x = [C12(1) D2(1) D1(1)  C11(1)];% 8 8 10 14 10 14 NaN 4 4 6 9 15];
y = [C12(2) D2(2) D1(2)  C11(2)];% 11 7  6 10 10  6 NaN 0 3 4 3  6];
%     mapshow(x,y,'Marker','+')
[xiA1,yiA1] = polyxpoly(x,y,xbox,ybox)
%     mapshow(xiA1,yiA1,'DisplayType','point','Marker','o')
x = [C12(1) B2(1) B1(1)  C11(1)];% 8 8 10 14 10 14 NaN 4 4 6 9 15];
y = [C12(2) B2(2) B1(2)  C11(2)];% 11 7  6 10 10  6 NaN 0 3 4 3  6];
%     mapshow(x,y,'Marker','+')
[xiA2,yiA2] = polyxpoly(x,y,xbox,ybox)
%     mapshow(xiA2,yiA2,'DisplayType','point','Marker','o')
%     hold off
%     hold on
%     % Fill 1
x = [C21(1) C1(1) C2(1)  C22(1)];% 8 8 10 14 10 14 NaN 4 4 6 9 15];
y = [C21(2) C1(2) C2(2)  C22(2)];% 11 7  6 10 10  6 NaN 0 3 4 3  6];
%     mapshow(x,y,'Marker','+')
[xiB1,yiB1] = polyxpoly(x,y,xbox,ybox);
%     mapshow(xiB1,yiB1,'DisplayType','point','Marker','o')
x = [C21(1) A1(1) A2(1)  C22(1)];% 8 8 10 14 10 14 NaN 4 4 6 9 15];
y = [C21(2) A1(2) A2(2)  C22(2)];% 11 7  6 10 10  6 NaN 0 3 4 3  6];
%     mapshow(x,y,'Marker','+')
[xiB2,yiB2] = polyxpoly(x,y,xbox,ybox);
%     mapshow(xiB2,yiB2,'DisplayType','point','Marker','o')
[xiA1 yiA1]= RescalePoly([xiA1 yiA1],Center,Factor);
[xiA2 yiA2]= RescalePoly([xiA2 yiA2],Center,Factor);
[xiB1 yiB1]= RescalePoly([xiB1 yiB1],Center,Factor);
[xiB2 yiB2]= RescalePoly([xiB2 yiB2],Center,Factor);

rectangle('Position',[Center(1)-L Center(2)-H 2*L 2*H],'FaceColor',[0 0 0]);
[Center_(:,1),Center_(:,2)]=RescalePoly([Center(1)-L Center(2)-H],Center,Factor);
rectangle('Position',[Center_(1) Center_(2) 2*L*(Factor) 2*H*(Factor)],'FaceColor',[1 1 1]);
fill([xiA1' xiA2'], [yiA1' yiA2'],rgb)
fill([xiA1' flip(xiA2')], [yiA1' flip(yiA2')],rgb)
fill([xiB1' flip(xiB2')], [yiB1' flip(yiB2')],rgb)
fill([xiB1' (xiB2')], [yiB1' (yiB2')],rgb)
% Fill 2
%     hold off
%     fill([B1,B2,D2,D1],'k')
end
