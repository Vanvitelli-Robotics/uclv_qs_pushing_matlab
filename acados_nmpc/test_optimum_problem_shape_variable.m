global plant S_p_x S_p_y

plant = p;
S_p_x = -0.04;
S_p_y = 0.025;

x0 = 0;
options = optimoptions('fminunc','Display','off','UseParallel',false);
% tic

% [x,fval,exitflag,output] = fminunc(@objfun,x0,options);
% toc
% disp(x)

%%
figure(1), hold on
ax1 = gca;
p.SP.printSpline(p.SP.a,p.SP.b,0.001);

figure(2), hold on
ax2 = gca;


for S_p_y=-0.04:0.01:0.04
    tic
    [x,fval,exitflag,output] = fminunc(@objfun,x0,options);
    toc

%     spl = p.SP.evalSpline(p.SP.FC,x);
%     plot(ax1,S_p_x,S_p_y,'*')
%     plot(ax1,spl(1),spl(2),'*')
%     figure(2), hold on
%     for j=-p.SP.b:0.001:p.SP.b
%         val = norm([S_p_x S_p_y] - plant.SP.evalSpline(plant.SP.FC,j))^2;
%         plot(ax2,j,val,'*b');
%     end
%     val = norm([S_p_x S_p_y] - plant.SP.evalSpline(plant.SP.FC,x0))^2;
%     plot(ax2,x0,val,'*r');
%     plot(ax2,x,fval,'*g');
%     figure(2), hold off
    x0 = x;
%     pause()
end
%%

% for S_p_y=-0.04:0.01:0.04
%     figure
%     for j=-p.SP.b:0.001:p.SP.b
%         val = norm([S_p_x S_p_y] - plant.SP.evalSpline(plant.SP.FC,j))^2;
%         plot(j,val,'*b');
%         hold on
%     end 
%     val = norm([S_p_x S_p_y] - plant.SP.evalSpline(plant.SP.FC,x0))^2;
%     plot(x0,val,'*r');
%     plot(x,fval,'*g');
%     grid on
%     x0 = x;
% end
% 
% figure
% p.SP.printSpline(p.SP.a,p.SP.b,0.001);
% hold on
% 
% spl = p.SP.evalSpline(p.SP.FC,x);
% plot(S_p_x,S_p_y,'*')
% plot(spl(1),spl(2),'*')

function f = objfun(x)
    global plant S_p_x S_p_y
        f = norm([S_p_x S_p_y] - plant.SP.evalSpline(plant.SP.FC,x))^2;
end