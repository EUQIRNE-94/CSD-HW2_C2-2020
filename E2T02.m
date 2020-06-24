% Ejercicio 2 - Tarea 2

close all
clear
clc

%% Sistema Transmisor de Rossler con método de identificcion de parametros lineal

tspan = [0 200];
dt = 0.1;

n = (tspan(2)/dt) +1;
tL = zeros(1,n);
lambda = zeros(1,n);
for i = 1:n
  tL(i) = (i-1)*dt;
  lambda(i) = LAMBDA(tL(i));
end

% Sistema de Rossler con cambio de variables
x0_R = [1,1,log(1),0,0,0,0,0,0,0,0,0,0,rand()];
[t,x] = ode45(@(t,x)lambdaRossler(t,x),tspan,x0_R);
x(:,3) = exp(x(:,3));

% Figuras Rossler con cambio de variables
figure(1)
set(gcf, 'Position', get(0, 'Screensize'));
subplot(3,1,1)
plot(t,x(:,1),'k','linewidth',2);
title('$x_1$','Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')

subplot(3,1,2)
plot(t,x(:,2),'k','linewidth',2);
title('$x_2$','Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')

subplot(3,1,3)
plot(t,x(:,3),'k','linewidth',2);
title({'$Senal - x_3$'},'Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')

set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E2_Estados.png')

% Retrato Fase
figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
plot3(x(:,1),x(:,2),x(:,3),'k','linewidth',1)
grid on
title('Retrato Fase Sistema','fontsize',30)
xlabel({'$x_1$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
zlabel({'$x_3$'},'Interpreter','latex','fontsize',20)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E2_RetratoFase.png')

% Lambda
figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
plot(tL,lambda,'g','linewidth',2); hold on; grid on
plot(t,x(:,13),'k','linewidth',2)
title('$\lambda \ y \ \hat{\lambda}$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$\lambda \ / \ \hat{\lambda}$'},'Interpreter','latex','fontsize',20)
legend({'$\lambda$','$\hat{\lambda}$'},'interpreter','latex','fontsize',16)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E2_Lambda.png')

% Error de Lambda
lambdaI = interp1(t,x(:,13),tL);
eL = lambda - lambdaI;
figure(4)
set(gcf, 'Position', get(0, 'Screensize'));
plot(tL,eL,'g','linewidth',2); grid on
title('Error de \lambda','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$error$'},'Interpreter','latex','fontsize',20)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E2_ErrorLambda.png')


%% Funciones
% Obtener lambda de la salida del sistema
function dx = lambdaRossler(t,x)
%% Sistema
x1 = x(1);
x2 = x(2);
x3 = x(3);
lambda = LAMBDA(t);

dx(1) = -x2 - exp(x3);
dx(2) = x1 + lambda*x2;
dx(3) = x1 + 2*exp(-x3) - 4;
y = x3;

%% Reconstruccion de Lambda
x3p = exp(x3);
u1 = -x3p;
u2  = (2/x3p) - 4;

k0 = 512;
k1 = 192;
k2 = 24;
v = 800;
gamma = 0.002;

w01 = x(4);
w02 = x(5);
w03 = x(6);
w11 = x(7);
w12 = x(8);
w13 = x(9);
w21 = x(10);
w22 = x(11);
w23 = x(12);
lambdag = x(13);
p = x(14);

dx(4) = w02;
dx(5) = w03;
dx(6) = -k0*w01 -k1*w02 - k2*w03 + y;
dx(7) = w12;
dx(8) = w13;
dx(9) = -k0*w11 -k1*w12 - k2*w13 + u1;
dx(10) = w22;
dx(11) = w23;
dx(12) = -k0*w21 -k1*w22 - k2*w23 + u2;

phi0 = k0*w01 + (k1 - 1)*w02 + k2*w03 + w12 + w21 + w23;
phi1 = w03 - w11 - w22;

yg = phi0 + lambdag*phi1;

dx(13) = -v*phi1*p*(yg - y);
dx(14) = -v*( (phi1^2)*(p^2) - gamma*p );

dx = dx';
end

function y = LAMBDA(t)
if t < 100
  y = 0.3;
else
  y = 0.3 + 0.2*sin((pi/25)*t);
end
end