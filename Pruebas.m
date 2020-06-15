clear
clc
close all

%% Sistema continuo
global sigma_L r_L b_L
sigma_L = 10;
r_L = 28;
b_L = 8/3;
tspan = [0 10];
x0_L = [1,1,1];
[t_L,x_L] = ode45(@LORENZ,tspan,x0_L);

%% Sistema Discreto
tf = 10;
dt = 0.01;
tspan = 0:dt:tf;
n = length(tspan);

sigma = 10;
r = 28;
b = 8/3;
A = [-sigma,sigma,0; r, -1, 0; 0,0,-b];
x(:,1) = rand(3,1);
for k = 1:n
  % Ruido del sistema y salida
  eps(:,k) = 0.1*randn(1,3);
  v(1,k) = 0.1*randn();
  % Calcular Sistema con ruido
  B = [0;-x(1,k)*x(3,k);x(1,k)*x(2,k)];
  x(:,k+1) = x(:,k) + ( A*x(:,k) + B )*dt + eps(:,k);
  y(1,k) = x(1,k) + v(1,k);
end

%% Graficas


%% Funciones
% Atractor de Lorenz
function dx = LORENZ(~,x)
global sigma_L r_L b_L

x1 = x(1);
x2 = x(2);
x3 = x(3);

dx(1) = sigma_L*(x2 - x1);
dx(2) = r_L*x1 - x2 - x1*x3;
dx(3) = -b_L*x3 + x1*x2;

dx = dx';
end