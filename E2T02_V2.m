% Ejercicio 2 - Tarea 2

close all
clear
clc

%% Sistema Transmisor de Rossler

tspan = [0 50];
dt = 0.1;

n = (tspan(2)/dt) +1;
tL = zeros(1,n);
lambda = zeros(1,n);
for i = 1:n
  tL(i) = (i-1)*dt;
  lambda(i) = LAMBDA(tL(i));
end

% Sistema de Rossler sin cambio de variables y observador parcial
x0_R = [1,1,1,0,0];
[t,x] = ode45(@(t,x)sisrossler(t,x),tspan,x0_R);

% Figuras Rossler sin cambio de variables
fig1(t,x,1)

% Sistema de Rossler con cambio de variables y observador parcial
x0_R = [1,1,log(1),0,0];
[tv,xv] = ode45(@(t,x)varrossler(t,x),tspan,x0_R);

% Figuras Rossler con cambio de variables
xv(:,3) = exp(xv(:,3));
fig1(tv,xv,4)

%% Funciones
% Obtener lambda de la salida del sistema
function dx = lambdaRossler(t,x)
%% Sistema
x1 = x(1);
x2 = x(2);
x3 = x(3);

lambda = LAMBDA(t);

dx(1) = -x2 - x3;
dx(2) = x1 + lambda*x2;
dx(3) = 2 + (x1 - 4)*x3;
y = log(x3);

%% Reconstruccion de Lambda
u0 = log(x3);
u1 = -x3;
u2  = (2/x3) - 4;



end

% Sistema con cambio de variables y observador parcial
function dx = varrossler(t,x)

x1 = x(1);
x2 = x(2);
x3 = x(3);
lambda = LAMBDA(t);

dx(1) = -x2 - exp(x3);
dx(2) = x1 + lambda*x2;
dx(3) = x1 + 2*exp(-x3) - 4;

% Observador parcial
xg1 = x(4);
xg2 = x(5);

dx(4) = -xg2 - exp(x3);
dx(5) = xg1 + lambda*xg2;

dx = dx';
end


function y = LAMBDA(t)
if t < 100
  y = 0.3;
else
  y = 0.3 + 0.2*sin((pi/25)*t);
end
end

% Sistema sin cambio de variables
function dx = sisrossler(t,x)
x1 = x(1);
x2 = x(2);
x3 = x(3);

lambda = LAMBDA(t);

dx(1) = -x2 - x3;
dx(2) = x1 + lambda*x2;
dx(3) = 2 + (x1 - 4)*x3;

% Observador parcial
xg1 = x(4);
xg2 = x(5);

dx(4) = -xg2 - x3;
dx(5) = -xg1 + lambda*xg2;

dx = dx';
end

% Figuras del sistema Rossler con observador parcial
function fig1(t,x,f)
figure(f)
subplot(3,1,1)
plot(t,x(:,1),'k','linewidth',2); hold on; grid on
plot(t,x(:,4),'g','linewidth',1)
title('$x_1 \ \ \hat{x}_1$','Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')
legend({'Sistema','Observador parcial'},'fontsize',16)

subplot(3,1,2)
plot(t,x(:,2),'k','linewidth',2); hold on; grid on
plot(t,x(:,5),'g','linewidth',1)
title('$x_2 \ \ \hat{x}_2$','Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')
legend({'Sistema','Observador parcial'},'fontsize',16)

subplot(3,1,3)
plot(t,x(:,3),'k','linewidth',2); hold on; grid on
plot(t,x(:,3),'g','linewidth',1)
title('$x_3 \ \ \hat{x}_3$','Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')
legend({'Sistema','Observador parcial'},'fontsize',16)

% Errores
e1 = x(:,1) - x(:,4);
e2 = x(:,2) - x(:,5);
e3 = x(:,3) - x(:,3);

figure(f+1)
plot(t,e1,'r',t,e2,'g',t,e3,'b','linewidth',2)
grid on
title('Errores','fontsize',30)
xlabel('Tiempo [t]','Interpreter','latex','fontsize',20)
ylabel({'error'},'Interpreter','latex','fontsize',20)
legend({'Error_1','Error_2','Error_3'},'fontsize',16)

% Retrato Fase
figure(f+2)
subplot(1,2,1)
plot3(x(:,1),x(:,2),x(:,3),'k','linewidth',3)
grid on
title('Retrato Fase Sistema','fontsize',30)
xlabel({'$x_1$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
zlabel({'$x_3$'},'Interpreter','latex','fontsize',20)

subplot(1,2,2)
plot3(x(:,4),x(:,5),x(:,3),'g','linewidth',1)
grid on
title('Retrato Fase Observador parcial','fontsize',30)
xlabel({'$x_1$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
zlabel({'$x_3$'},'Interpreter','latex','fontsize',20)
end