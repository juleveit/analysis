function plot_orrf_deg(gf,nsigma,style,width)

%TAKE CARE: angle theta is measured from x being vertical to compare to
%grating where 0deg = vertical
%this is why x and y are exchanged for the spread. for center they are
%correct!

% rotation matrix to rotate the axes with respect to an angle theta
R = [ cos(gf.theta) sin(gf.theta); -sin(gf.theta) cos(gf.theta) ];


%rotate in origin, then shift to correct center of rf
ver_line = [ [0 0]; nsigma*gf.yspreadDeg*[-1,1] ];
horz_line = [ nsigma*gf.xspreadDeg*[-1 1]; [0, 0] ];
new_ver_line    = R*ver_line;
new_horz_line   = R*horz_line;
new_ver_line = new_ver_line + [gf.xcenterDeg, gf.xcenterDeg ; gf.ycenterDeg, gf.ycenterDeg];
new_horz_line = new_horz_line + [gf.xcenterDeg, gf.xcenterDeg ; gf.ycenterDeg, gf.ycenterDeg];

% the ellipse
%rotate in origin, then shit to rf center
theta_r = linspace(0,2*pi);
ellipse_x = nsigma*gf.xspreadDeg*cos(theta_r);
ellipse_y = nsigma*gf.yspreadDeg*sin(theta_r);
rotated_ellipse = R * [ellipse_x;ellipse_y];
rotated_ellipse(1,:) = rotated_ellipse(1,:)+gf.xcenterDeg;
rotated_ellipse(2,:) = rotated_ellipse(2,:)+gf.ycenterDeg;

%plot
hold on
plot( new_ver_line(1,:),new_ver_line(2,:),style,'LineWidth',width );
plot( new_horz_line(1,:),new_horz_line(2,:),style,'LineWidth',width );
plot( rotated_ellipse(1,:),rotated_ellipse(2,:),style,'LineWidth',width );