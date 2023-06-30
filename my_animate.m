
function my_animate(x, y, theta)
% Function to animate the trajectory tracking of the pusher slider system
% Input: 2D position [x,y,theta]
    p0 = [x(1) y(1) 0];
    quat0 = quaternion(rotz(theta(1)),'rotmat','frame');
    patch = poseplot(quat0,p0);
    patch.ScaleFactor = 0.1;
    ax1 = gca; ax1.View = [0.7528 90];
    xlim(ax1,[-0.1 0.5])
    ylim(ax1,[-0.3 0.3])
    xlabel("x [m]")
    ylabel("y [m]")
    zlabel("z [m]")

    for i = 2:1:length(x)
        position = [x(i) y(i) 0];
        quat = quaternion(rotz(rad2deg(theta(i))),'rotmat','frame');
        set(patch,Orientation=quat,Position=position);
        pause(0.1)
        drawnow
    end
