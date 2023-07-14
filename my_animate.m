
function my_animate(x, y, theta, rx, ry,sampling_time, traj)
    % Function to animate the trajectory tracking of the pusher slider system
    % Input: 2D position [x,y,theta]
    
    % Slider frame
    %     theta = -theta;
    p0_s = [x(1) y(1) 0];
    quat0_s = quaternion(helper.my_rotz(theta(1)),'rotmat','frame');
    slider_p = poseplot(quat0_s,p0_s, MeshFileName="cad_models/cuboide_santal.stl", ScaleFactor=0.001);%,PatchFaceColor="yellow");
    hold on
    
    % Pusher frame
    %     p0_p = helper.my_rotz(theta(1))*[rx(1) ry(1) 0]';
    %     quat0_p = quaternion(helper.my_rotz(0),'rotmat','frame');
    %     pusher_p = poseplot(quat0_p,p0_p, MeshFileName="cad model pusher\cad_model_pusher_1.stl",ScaleFactor=0.001);
    radius_pusher = 0.003;
    pusher_position = [x(1) y(1) 0]' + helper.my_rotz(theta(1))*[rx(1) ry(1) 0]';
    viscircles([pusher_position(1) pusher_position(2)],radius_pusher,"Color",'red',"LineWidth",6);
    
    hold off
    
    ax1 = gca; ax1.View = [-0.1303  -90];
    xlim(ax1,[-0.1 0.5])
    ylim(ax1,[-0.3 0.3])
    xlabel("x [m]")
    ylabel("y [m]")
    zlabel("z [m]")
    
    for i = 2:1:length(x)
        position_s = [x(i) y(i) 0]';
        quat_s = quaternion([0 0 theta(i)],'rotvec');
        set(slider_p,Orientation=quat_s,Position=position_s);
        hold on
        pusher_position = [x(i) y(i) 0]'+ helper.my_rotz(theta(i))*[rx(i) ry(i) 0]';
        viscircles([pusher_position(1) pusher_position(2)],radius_pusher,"Color",'red',"LineWidth",6);
        plot(traj(1,:), traj(2,:), '-.','LineWidth',3,'Color','b');
        pause(sampling_time)
        drawnow
    end
end
