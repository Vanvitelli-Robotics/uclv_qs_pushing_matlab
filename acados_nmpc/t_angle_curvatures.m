s_values = p.SP.a:0.001:p.SP.b;
t_vers = zeros(length(s_values),2);

t_vers_ = p.SP.FC_dot(1);
% t_vers = full(t_vers_(s_values));

for i = 1:length(s_values)
    t_vers(i,:) = full(t_vers_(s_values(i)));
end

t_diff = diff(t_vers./vecnorm(t_vers))./diff(s_values');
t_diff = vecnorm(t_diff');

t_norm = vecnorm(t_vers');

t_angle = atan2(t_vers(:,2),t_vers(:,1));


t_curvatures = diff(t_angle');%./diff(s_values);

t_curvatures(abs(t_curvatures) > 3*pi/2) = t_curvatures(abs(t_curvatures) > 3*pi/2) - 2*pi*sign(t_curvatures(abs(t_curvatures) > 3*pi/2));

t_curvatures = t_curvatures./diff(s_values);

curvatures = p.SP.getCurvatures(s_values);

figure, plot(s_values(1:end-1), abs(t_curvatures/max(t_curvatures))), hold on
figure, plot(s_values(1:end-1), t_diff/max(t_diff))
figure, plot(s_values, curvatures/max(curvatures))


%%
import casadi.*

s_values = p.SP.a:0.0001:p.SP.b;

s = SX.sym('s');

t_vers_ = p.SP.FC_dot(1);
t_vers_s = t_vers_(s);

t_angle = atan2(t_vers_s(:,2),t_vers_s(:,1));


t_angle_dot = Function('t_angle_dot',{s},{gradient(t_angle,s)});

t_angle_dot_values = [];
for j = s_values
    t_angle_dot_values = [t_angle_dot_values full(t_angle_dot(j))];
end






t_curvatures = diff(t_angle');%./diff(s_values);

t_curvatures(abs(t_curvatures) > 3*pi/2) = t_curvatures(abs(t_curvatures) > 3*pi/2) - 2*pi*sign(t_curvatures(abs(t_curvatures) > 3*pi/2));

t_curvatures = t_curvatures./diff(s_values);

%%

s_values = p.SP.a:0.0001:p.SP.b;
v_boundd = zeros(1,length(s_values));
t_angle = zeros(1,length(s_values));
for i = 1 : length(s_values)
    [v_boundd(i), t_angle(i)] = controller.update_tangential_velocity_bounds(s_values(i));
end

% obj_txt = 'montana';
% 
% open('test_curvatures/objects.fig')
subplot(3,1,1),hold on, plot(s_values,v_boundd)%;,'DisplayName',obj_txt)
subplot(3,1,2),hold on, plot(s_values,-t_angle);%,'DisplayName',obj_txt)
t_angle = 0:0.01:250;
subplot(3,1,3),hold on, plot(t_angle, min(0.5*controller.v_alpha./(abs(t_angle-3)+0.0001)+controller.d_v_bound, controller.u_t_ub));%,'DisplayName',obj_txt)

%% 
s_values = params.x_sim(4,:);
v_boundd = zeros(1,length(s_values));
t_angle = zeros(1,length(s_values));
for i = 1 : length(s_values)
    [v_boundd(i), t_angle(i)] = controller.update_tangential_velocity_bounds(s_values(i));
end
figure
subplot(3,1,1),hold on, plot(params.t,v_boundd)
subplot(3,1,2),hold on, plot(params.t,t_angle)
