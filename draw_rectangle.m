function[]= draw_rectangle(center_location,L,H,theta,rgb,Num)
center1=center_location(1);
center2=center_location(2);
R= ([cos(theta), -sin(theta), sin(theta), cos(theta)]);
X=([-L/2, L/2, L/2, -L/2]);
Y=([-H/2, -H/2, H/2, H/2]);

x_lower_left=center1+ (X(1)*cos(theta)-Y(1)*sin(theta));
x_lower_right=center1+(X(2)*cos(theta)-Y(2)*sin(theta));
x_upper_right=center1+(X(3)*cos(theta)-Y(3)*sin(theta));
x_upper_left=center1+(X(4)*cos(theta)-Y(4)*sin(theta));
y_lower_left=center2+(X(1)*sin(theta)+Y(1)*cos(theta));
y_lower_right=center2+(X(2)*sin(theta)+Y(2)*cos(theta));
y_upper_right=center2+(X(3)*sin(theta)+Y(3)*cos(theta));
y_upper_left=center2+(X(4)*sin(theta)+Y(4)*cos(theta));
x_coor=[x_lower_left x_lower_right x_upper_right x_upper_left];
y_coor=[y_lower_left y_lower_right y_upper_right y_upper_left];
fill(x_coor, y_coor,rgb,'EdgeColor','white','LineWidth',Num);
axis equal;
end
