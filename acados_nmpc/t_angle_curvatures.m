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
plot(s_values(1:end-1), t_diff/max(t_diff))
plot(s_values, curvatures/max(curvatures))


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
curv = zeros(1,length(s_values));
t_angle = zeros(1,length(s_values));
for i = 1 : length(s_values)
    [curv(i), t_angle(i)] = controller.update_tangential_velocity_bounds(s_values(i));
end

figure, plot(s_values,curv)
figure, plot(s_values,t_angle)

