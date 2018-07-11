function plot_orrf_pixcorr(gf,nsigma,style,width)

    %rotation matrix other way around because with rf plot axis is ij

    % rotation matrix to rotate the axes with respect to an angle theta
    R = [ cos(gf.theta) sin(gf.theta); -sin(gf.theta) cos(gf.theta) ];

    if gf.yspreadDeg == gf.shortedge
        help = gf.yspreadPix;
        gf.yspreadPix = gf.xspreadPix;
        gf.xspreadPix = help;
    end

    %rotate in origin, then shift to correct center of rf
    ver_line = [ [0 0]; nsigma*gf.xspreadPix*[-1,1] ];
    horz_line = [ nsigma*gf.yspreadPix*[-1 1]; [0, 0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    new_ver_line = new_ver_line + [gf.xcenterPix, gf.xcenterPix ; gf.ycenterPix, gf.ycenterPix];
    new_horz_line = new_horz_line + [gf.xcenterPix, gf.xcenterPix ; gf.ycenterPix, gf.ycenterPix];

    % the ellipse
    %rotate in origin, then shit to rf center
    theta_r = linspace(0,2*pi);
    ellipse_x = nsigma*gf.xspreadPix*cos(theta_r);
    ellipse_y = nsigma*gf.yspreadPix*sin(theta_r);
    rotated_ellipse = R * [ellipse_y;ellipse_x];
    rotated_ellipse(1,:) = rotated_ellipse(1,:)+gf.xcenterPix;
    rotated_ellipse(2,:) = rotated_ellipse(2,:)+gf.ycenterPix;

    %plot
    hold on
    plot( new_ver_line(1,:),new_ver_line(2,:),style,'LineWidth',width );
    plot( new_horz_line(1,:),new_horz_line(2,:),style,'LineWidth',width );
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),style,'LineWidth',width );
end