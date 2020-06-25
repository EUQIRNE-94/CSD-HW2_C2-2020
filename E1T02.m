% Ejercicio 1 - Tarea 2
close all
clear
clc

%% Sistema de Lorenz sujeto a ruido blanco en proceso y medición

% Tiempo de simulación
tf = 5;
dt = 0.01;
tspan = 0:dt:tf;
n = length(tspan);

%% Método Discreto - FKE general

% Prealocacion de variables
x = zeros(3,n);                                                             % Sistema
xg = zeros(3,n);                                                            % Observador del sistema
y = zeros(1,n);                                                             % Salida del sistema
xg_k = zeros(3,n);                                                          % x(k+1|k)
K = zeros(3,n);                                                             % Ganancias
P = zeros(3,3,n);                                                           % Error de covarianza
Q = zeros(3,3,n);                                                           % Matriz de covarianza del sistema
R = zeros(1,n);                                                             % Covarianza de la salida
eps = zeros(3,n);                                                           % Ruidos del proceso
v = zeros(1,n);                                                             % Ruido de lectura
e1 = zeros(1,n);                                                            % Error de variable medida

% Ruido del proceso y lectura
mu = 0;                                                                     % Centro del ruido
sig = 0.01;                                                                 % Varianza
% Condiciones iniciales
x(:,1) = rand(3,1);
xg_k(:,1) = x(:,1);
P(:,:,1) = rand(3);

% Parametros del sistema
sigma = 10;
r = 28;
b = 8/3;

% Matrices lineales
A = [-sigma,sigma,0; r, -1, 0; 0,0,-b];

for k = 1:n
  % Ruido del sistema y salida
  eps(:,k) = mu + sig*randn(1,3);
  v(1,k) = mu + sig*randn(1,1);
  % Calcular Sistema con ruido
  B = [0;-x(1,k)*x(3,k);x(1,k)*x(2,k)];
  x(:,k+1) = x(:,k) + ( A*x(:,k) + B )*dt + eps(:,k);
  y(1,k) = x(1,k) + v(1,k);
  
  % Obtener el estimado a posteriori del estado
  Bg = [0;-xg_k(1,k)*xg_k(3,k);xg_k(1,k)*xg_k(2,k)];
  xg_k(:,k+1) = xg_k(:,k) + ( A*xg_k(:,k) + Bg )*dt;
  
  % Obtener error de covarianza estimado
  F = [ 1-sigma*dt,sigma*dt,0 ; (r-x(3,k))*dt,1-dt,-x(1,k)*dt ; x(2,k)*dt,x(1,k)*dt,1-b*dt ];
  c_e1 = covarianza(eps(1,1:k),eps(1,1:k));
  c_e2 = covarianza(eps(2,1:k),eps(2,1:k));
  c_e3 = covarianza(eps(3,1:k),eps(3,1:k));
%   Q(:,:,k) = [c_e1,0,0;0,c_e2,0;0,0,c_e3];
  Q(:,:,k) = [sig*sig,0,0;0,sig*sig,0;0,0,sig*sig];
  P_k = F*P(:,:,k)*F' + Q(:,:,k);
  
  % Obtener ganancias de Kalman
  H = [1,0,0];
%   R(1,k) = covarianza(v(1,1:k),v(1,1:k));
  R(1,k) = sig*sig;
  K(:,k) = P_k*H'*inv( H*P_k*H' + R(1,k) );
  
  % Generar copia del sistema con ganancias de Kalman
  e1(1,k) = y(1,k) - xg(1,k);
  xg(:,k+1) = xg_k(:,k) + K(:,k)*(e1(1,k))*dt;
  
  % Actualizar error de covarianza
  P(:,:,k+1) = (eye(3) - K(:,k)*H)*P_k;  
end

% Errores
e2 = x(2,:) - xg(2,:);
e3 = x(3,:) - xg(3,:);

%% Metodo continuo - Paper
tspan_c = [0 tf];

x0_L = [x(1,1),x(2,1),x(3,1),0,0,0,.1,.1,.1,.1,.1,.1,.1,.1,.1];
[t_L,x_L] = ode45(@LORENZ,tspan_c,x0_L);

%% Figuras
figure(1)
set(gcf, 'Position', get(0, 'Screensize'));
% Estado 1
subplot(3,1,1)
plot(tspan,x(1,1:n),'k','linewidth',4);hold on;grid on
plot(tspan,xg(1,1:n),'g','linewidth',2)
title('Estado $x_1$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_1$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% Estado 2
subplot(3,1,2)
plot(tspan,x(2,1:n),'k','linewidth',4);hold on; grid on
plot(tspan,xg(2,1:n),'g','linewidth',2)
title('Estado $x_2$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% Estado 3
subplot(3,1,3)
plot(tspan,x(3,1:n),'k','linewidth',4);hold on; grid on
plot(tspan,xg(3,1:n),'g','linewidth',2)
title('Estado $x_3$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_3$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'E1_Estados_Disc.png')

% Errores
figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
plot(tspan,e1,'r','linewidth',2); hold on; grid on
plot(tspan,e2(1:n),'g','linewidth',2)
plot(tspan,e3(1:n),'b','linewidth',2)
title('Errores','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$e_i$'},'Interpreter','latex','fontsize',20)
legend({'e_1','e_2','e_3'},'fontsize',16)

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'E1_Errores_Disc.png')

% Retrato Fase
figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
plot3(x(1,:),x(2,:),x(3,:),'k','linewidth',1); hold on; grid on
plot3(xg(1,:),xg(2,:),xg(3,:),'g','linewidth',1)
title('Retrato Fase','fontsize',30)
xlabel({'$x_1$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
zlabel({'$x_3$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'E1_RetratoFase_Disc.png')

%% Figuras sistema continuo
figure(4)
set(gcf, 'Position', get(0, 'Screensize'));
% Estado 1
subplot(3,1,1)
plot(t_L,x_L(:,1),'k','linewidth',4);hold on;grid on
plot(t_L,x_L(:,4),'g','linewidth',2);
title('Estado $x_1$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_1$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% Estado 2
subplot(3,1,2)
plot(t_L,x_L(:,2),'k','linewidth',4);hold on;grid on
plot(t_L,x_L(:,5),'g','linewidth',2);
title('Estado $x_2$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% Estado 3
subplot(3,1,3)
plot(t_L,x_L(:,3),'k','linewidth',4);hold on;grid on
plot(t_L,x_L(:,6),'g','linewidth',2);
title('Estado $x_3$','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_3$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'E1_Estados_Cont.png')

% Errores
e1L = x_L(:,1) - x_L(:,4);
e2L = x_L(:,2) - x_L(:,5);
e3L = x_L(:,3) - x_L(:,6);

figure(5)
set(gcf, 'Position', get(0, 'Screensize'));
plot(t_L,e1L,'r','linewidth',2); hold on; grid on
plot(t_L,e2L,'g','linewidth',2)
plot(t_L,e3L,'b','linewidth',2)
title('Errores','interpreter','latex','fontsize',30)
xlabel({'Tiempo $t$'},'Interpreter','latex','fontsize',20)
ylabel({'$e_i$'},'Interpreter','latex','fontsize',20)
legend({'e_1','e_2','e_3'},'fontsize',16)

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'E1_Errores_Cont.png')

% Retrato Fase
figure(6)
set(gcf, 'Position', get(0, 'Screensize'));
plot3(x_L(:,1),x_L(:,2),x_L(:,3),'b','linewidth',1); hold on; grid on
plot3(x_L(:,4),x_L(:,5),x_L(:,6),'g','linewidth',1);
title('Retrato Fase','fontsize',30)
xlabel({'$x_1$'},'Interpreter','latex','fontsize',20)
ylabel({'$x_2$'},'Interpreter','latex','fontsize',20)
zlabel({'$x_3$'},'Interpreter','latex','fontsize',20)
legend({'Sistema','Observador',},'fontsize',16)

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'E1_RetratoFase_Cont.png')

%% Funciones
% FKE
function c = covarianza(x,y)
n_x = length(x);
n_y = length(y);

if n_x == n_y
  % Promedio de variables
  x_b = mean(x);
  y_b = mean(y);
  % Calculo de covarianza
  c = sum( (x - x_b).*(y - y_b)  )/n_x;
end

end

% Paper
function dx = LORENZ(~,x)

sigma = 10;
r = 28;
b = 8/3;
sig = 0.01;

x1 = x(1);
x2 = x(2);
x3 = x(3);

E1 = sig*randn();
E2 = sig*randn();
E3 = sig*randn();

dx(1) = sigma*(x2 - x1) + E1;
dx(2) = r*x1 - x2 - x1*x3 + E2;
dx(3) = -b*x3 + x1*x2 + E3;

v = sig*randn();
y = x1 + v;

% Filtro Kalman
p11 = x(7);
p12 = x(8);
p13 = x(9);
p21 = x(10);
p22 = x(11);
p23 = x(12);
p31 = x(13);
p32 = x(14);
p33 = x(15);

P = [p11,p12,p13;p21,p22,p23;p31,p32,p33];
F = [-sigma,sigma,0 ; r-x3,1,-x1 ; x2,x1,-b];
H = [1,0,0];
Q = [sig^2,0,0;0,sig^2,0;0,0,sig^2];
R = sig^2;

Pp = F*P + P*F' - P*H'*(1/R)*H*P + Q;

dx(7) = Pp(1,1);
dx(8) = Pp(1,2);
dx(9) = Pp(1,3);
dx(10) = Pp(2,1);
dx(11) = Pp(2,2);
dx(12) = Pp(2,3);
dx(13) = Pp(3,1);
dx(14) = Pp(3,2);
dx(15) = Pp(3,3);

K = P*[1;0;0]*(1/R);

xg1 = x(4);
xg2 = x(5);
xg3 = x(6);

e1 = y - xg1;

dx(4) = sigma*(xg2 - xg1) + K(1)*e1;
dx(5) = r*xg1 - xg2 - xg1*xg3 + K(2)*e1;
dx(6) = -b*xg3 + xg1*xg2 + K(3)*e1;

dx = dx';
end