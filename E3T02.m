% Ejercicio 3 - Tarea 2
close all
clear
clc

%% Sistema Transmisor de Rossler con observador adaptable para comunicaciones
tspan = [0 2];
dt = 0.01;
tspan1 = 0:dt:2;

n = (tspan(2)/dt) +1;
tL = zeros(1,n);
lambda = zeros(1,n);
for i = 1:n
  tL(i) = (i-1)*dt;
  lambda(i) = LAMBDA(tL(i));
end

% Sistema de Rossler con cambio de variables
x0_R = [.1,.1,log(.1),.1,.1,log(.1),.1,.1,0,0,0,0,0,0,0.1];
[t,x] = ode45(@(t,x)lambdaRossler(t,x),tspan1,x0_R);
x(:,3) = exp(x(:,3));

% Figuras Rossler con cambio de variables
% figure(1)
% set(gcf, 'Position', get(0, 'Screensize'));
% subplot(3,1,1)
% plot(t,x(:,1),'k','linewidth',2);hold on
% plot(t,x(:,4),'r','linewidth',2);
% plot(t,x(:,7),'b','linewidth',2);
% title('$x_1$','Interpreter','latex','fontsize',30)
% xlabel('Tiempo [t]')
% ylabel('Uds')
% legend({'$\xi$','Z','$\xi_{filtrada}$'},'interpreter','latex','Fontsize',16)
% 
% subplot(3,1,2)
% plot(t,x(:,2),'k','linewidth',2);hold on
% plot(t,x(:,5),'r','linewidth',2);
% plot(t,x(:,8),'b','linewidth',2);
% title('$x_2$','Interpreter','latex','fontsize',30)
% xlabel('Tiempo [t]')
% ylabel('Uds')
% legend({'$\xi$','Z','$\xi_{filtrada}$'},'interpreter','latex','Fontsize',16)
% 
% subplot(3,1,3)
% plot(t,x(:,3),'k','linewidth',2);hold on
% plot(t,x(:,6),'r','linewidth',2);
% title({'$Senal - x_3$'},'Interpreter','latex','fontsize',30)
% xlabel('Tiempo [t]')
% ylabel('Uds')
% legend({'$\xi$','Z'},'interpreter','latex','Fontsize',16)
% 
% % Retrato Fase
% figure(2)
% set(gcf, 'Position', get(0, 'Screensize'));
% plot3(x(:,1),x(:,2),x(:,3),'k','linewidth',1);hold on
% plot3(x(:,4),x(:,5),x(:,6),'r','linewidth',1)
% grid on
% title('Retrato Fase Sistema X y Z','fontsize',30)
% xlabel({'$x_1$'},'Interpreter','latex','fontsize',20)
% ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
% zlabel({'$x_3$'},'Interpreter','latex','fontsize',20)

figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
subplot(3,1,1)
plot(t,x(:,9),'k','linewidth',2);hold on
plot(t,x(:,12),'g','linewidth',1)
title({'$\eta_1 \ / \ \hat{\eta}_1$'},'Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')
legend({'$\eta_1$','$\hat{\eta}_1$'},'interpreter','latex','Fontsize',16)

subplot(3,1,2)
plot(t,x(:,10),'k','linewidth',2);hold on
plot(t,x(:,13),'g','linewidth',1);
title({'$\eta_1 \ / \ \hat{\eta}_1$'},'Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')
legend({'$\eta_1$','$\hat{\eta}_1$'},'interpreter','latex','Fontsize',16)


subplot(3,1,3)
plot(t,x(:,11),'k','linewidth',2);hold on
plot(t,x(:,14),'g','linewidth',1);
title({'$y \ / \ \hat{y}$'},'Interpreter','latex','fontsize',30)
xlabel('Tiempo [t]')
ylabel('Uds')
legend({'$y$','$\hat{y}$'},'interpreter','latex','Fontsize',16)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E3_Estados.png')

% Lambda
figure(4)
set(gcf, 'Position', get(0, 'Screensize'));
plot(tL,lambda,'g','linewidth',2); hold on; grid on
plot(t,x(:,15),'k','linewidth',2)
title('$\lambda \ y \ \hat{\lambda}$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$\lambda \ / \ \hat{\lambda}$'},'Interpreter','latex','fontsize',20)
legend({'$\lambda$','$\hat{\lambda}$'},'interpreter','latex','fontsize',16)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E3_Lambda.png')

% Error de Lambda
set(gcf, 'Position', get(0, 'Screensize'));
lambdaI = interp1(t,x(:,15),tL);
eL = lambda - lambdaI;
figure(5)
set(gcf, 'Position', get(0, 'Screensize'));
plot(tL,eL,'g','linewidth',2); grid on
title('Error de \lambda','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$error$'},'Interpreter','latex','fontsize',20)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'E3_ErrorLambda.png')

%% Funciones
% Obtener lambda de la salida del sistema
function dx = lambdaRossler(t,x)
%% Sistema
E1 = x(1);
E2 = x(2);
E3 = x(3);
lambda = LAMBDA(t);

% Sistema con cambio de variable
dx(1) = -E2 - exp(E3);
dx(2) = E1 + lambda*E2;
dx(3) = E1 + 2*exp(-E3) - 4;
y = E3;

% Cambio de coordenadas a Z
z1 = x(4);
z2 = x(5);
z3 = x(6);

dx(4) = 4*exp(-y) - 2 + lambda*exp(y);
dx(5) = z1 - z3 - exp(y) + lambda*(2 - 4*exp(-y));
dx(6) = z2 + 4*exp(-y) - 2 + lambda*y;
yz = z3;

% Transformacion Filtrada
k1 = 1;
k2 = 1;
E1f = x(7);
E2f = x(8);

dx(7) = k1*E2f + k1*yz + exp(yz);
dx(8) = E1f + k2*E2f + k2*yz - 4*exp(yz) + 2;

% Sistema con nuevas coordenadas
eta1 = x(9);
eta2 = x(10);
y_nc = x(11); %

dx(9) = k1*eta2 - k1*k2*y_nc + (k1 + 1)*(4*exp(-y_nc) - 2);
dx(10) = eta1 + k2*eta2 - ( k1 + k2*k2 + 1 )*y_nc + k2*(4*exp(-y_nc) - 2) - exp(y_nc);
dx(11) = eta2 + k2*y_nc + 4*exp(-y_nc) - 2 + lambda*(E2f + y_nc);  %

% Observador adaptable
l1 = .01;
l2 = .01;
l3 = .01;
gamma = 1;
eta1_g = x(12);
eta2_g = x(13);
y_ncg = x(14);
teta = x(15);

dx(12) = k1*eta2_g - k1*k2*y_ncg + (k1 + 1)*(4*exp(-y_nc) - 2) + l1*(y_ncg - y_nc);
dx(13) = eta1_g + k2*eta2_g - ( k1 + k2*k2 + 1 )*y_ncg + k2*(4*exp(-y_nc) - 2) - exp(y_nc) + l2*(y_ncg - y_nc);
dx(14) = eta2_g + k2*y_ncg + 4*exp(-y_nc) - 2 + teta*(E2f + y_nc) + l3*(y_ncg - y_nc);

dx(15) = -gamma*(E2f + y_nc)*(y_ncg - y_nc);

dx = dx';
end

function y = LAMBDA(t)
if t < 25
  y = 0.3;
else
  y = 0.3 + 0.2*sin((pi/25)*t - 25);
end
end