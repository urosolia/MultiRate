clc
clear all
close all

dtStr = '5';
cnstr = 0;
maxTime = 4;

dtStr = '1';
cnstr = 1;
maxTime = 4;

sample = 1;

if cnstr == 1;
    cnstrString = '_withConstr.mat';
else
    cnstrString = '.mat';
end;

load(strcat('dt_mpc',dtStr,'_withoutBarrier',cnstrString))
x_cl_withoutBarrier = x_cl;
xn_cl_withoutBarrier = xn_cl;
u_cl_withoutBarrier = u;
h_withoutBarrier = h;
V_withoutBarrier = V;
time_withoutBarrier = time;

%%
load(strcat('dt_mpc',dtStr,'_withBarrier',cnstrString))
x_cl_withBarrier = x_cl;
xn_cl_withBarrier = xn_cl;
u_cl_withBarrier = u;
h_withBarrier = h;
V_withBarrier = V;
time_withBarrier = time;
T = maxTime;

%% Plot states
figure 
subplot(4,1,1)
plot([time_withBarrier(1), time_withBarrier(end)],[0,0],'-k','linewidth',1);
ylim([-1.25,0.25])
hold on
plot(time_withoutBarrier, x_cl_withoutBarrier(1,:),'-b','linewidth',2);
plot(time_withBarrier, x_cl_withBarrier(1,:),'-r','linewidth',2);
xlim([0 T])
ylabel('x')
ylabel('$$x$$ [m]', 'interpreter','latex', 'fontsize',20)
subplot(4,1,2)
plot(time_withoutBarrier, x_cl_withoutBarrier(2,:),'-b','linewidth',2);
hold on
plot([time_withBarrier(1), time_withBarrier(end)],[0,0],'-k','linewidth',1);
plot(time_withBarrier, x_cl_withBarrier(2,:),'-r','linewidth',2);
xlim([0 T])
ylabel('$$v$$ [m/s]', 'interpreter','latex', 'fontsize',20)
subplot(4,1,3)
plot(time_withoutBarrier, x_cl_withoutBarrier(3,:),'-b','linewidth',2);
hold on
plot([time_withBarrier(1), time_withBarrier(end)],[0,0],'-k','linewidth',1);
plot(time_withBarrier, x_cl_withBarrier(3,:),'-r','linewidth',2);
xlim([0 T])
ylabel('$$\phi$$ [rad]', 'interpreter','latex', 'fontsize',20)
subplot(4,1,4)
b = plot(time_withoutBarrier, x_cl_withoutBarrier(4,:),'-b','linewidth',2);
hold on
d = plot([time_withBarrier(1), time_withBarrier(end)],[0,0],'-k','linewidth',1);
c = plot(time_withBarrier, x_cl_withBarrier(4,:),'-r','linewidth',2);
xlim([0 T])
h = legend([b,c, d], {'Linear MPC', 'Linear MPC with CLF-CBF', 'Goal state'});
set(h, 'interpreter','latex','fontsize',18)
ylabel('$$\dot \theta$$ [rad/s]', 'interpreter','latex', 'fontsize',20)
xlabel('time [s]', 'interpreter','latex', 'fontsize',20)



%%
figure
hold on
b = plot(time_withoutBarrier, x_cl_withoutBarrier(1,:),'-b','linewidth',2);
c = plot(time_withBarrier, x_cl_withBarrier(1,:),'-r','linewidth',2);
a = plot([time_withBarrier(1), time_withBarrier(end)],[0,0],'-k','linewidth',2);
xlim([0 T])
ylabel('x')
h = legend([a,b,c], {'Goal State', 'Linear MPC', 'Linear MPC with CLF-CBF'});
set(h, 'interpreter','latex','fontsize',18)
ylabel('$$x$$[m]', 'interpreter','latex', 'fontsize',20)
xlabel('time [s]', 'interpreter','latex', 'fontsize',20)

%%
if cnstr == 1
    figure
    hold on
    b = plot(time_withoutBarrier(1:sample:end), x_cl_withoutBarrier(3,1:sample:end),'-b','linewidth',2);
    c = plot(time_withBarrier(1:sample:end), x_cl_withBarrier(3,1:sample:end),'-r','linewidth',2);
    a = plot([time_withBarrier(1), time_withBarrier(end)],[0.78,0.78],'-k','linewidth',2);
    xlim([0 T]);
    h = legend([a,b,c], {'Constraint', 'Linear MPC', 'Linear MPC with CLF-CBF'});
    set(h, 'interpreter','latex','fontsize',18)
    ylabel('$$\theta$$ [rad]', 'interpreter','latex', 'fontsize',20)
    xlabel('time [s]', 'interpreter','latex', 'fontsize',20)
end
%% 
for i = 1:(round(1000/(2*pi))+1)
    xEl(i) = 0.1*sin(2*pi/round(1000/(2*pi))*i);
    yEl(i) = 0.1*cos(2*pi/round(1000/(2*pi))*i);
   
end
figure
a = plot(xEl, yEl,'-k','linewidth',2);
hold on
err_withoutBarrier = xn_cl_withoutBarrier - x_cl_withoutBarrier;
err_withBarrier = xn_cl_withBarrier - x_cl_withBarrier;
b = plot(err_withoutBarrier(2,:), err_withoutBarrier(4,:),'-b','linewidth',2);
c = plot(err_withBarrier(2,:), err_withBarrier(4,:),'-r','linewidth',2);
h = legend([a,b,c], {'$\delta S_e$', 'Linear MPC', 'Linear MPC with CLF-CBF'});
set(h, 'interpreter','latex','fontsize',18)
ylabel('$$\omega - \bar \omega$$ [rad]', 'interpreter','latex', 'fontsize',20)
xlabel('$v - \bar v$ [m/s]', 'interpreter','latex', 'fontsize',20)
ylim([-0.4, 0.25])
xlim([-0.11, 0.11])
%%
load(strcat('dt_mpc',dtStr,'_withBarrier',cnstrString))
T = maxTime;

figure
hold on
a = plot(time(1:end-1), uBarrier,'-b','linewidth',2);
b = plot(time(1:end-1), uMPC,'-r','linewidth',2);
c = plot(time(1:end-1), u,'-k','linewidth',2);
h = legend([a,b,c], {'$u$', '$v$', 'Total input $u+v$'});
set(h, 'interpreter','latex','fontsize',18)
ylabel('voltage [V]', 'interpreter','latex', 'fontsize',20)
xlabel('time [s]', 'interpreter','latex', 'fontsize',20)
xlim([0 T]);

