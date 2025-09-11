clear; close all; clc;
y(1:59) = 0; 
u(1:59) = 0;

%% Referência
yr = [0*ones(1,240)  2*ones(1,240)  4*ones(1,240)  2*ones(1,240)  1*ones(1,240)];

%% Parâmetros do controle
N  = 40;     % horizonte de predição
N1 = 1;      % início de predição
Ny = 10;     % horizonte de saída 
Nu = 6;      % horizonte de controle
lambda = 1;  % penalização do esforço de controle

%% Planta
nump = [1]; 
denp = [10 1];
Gp   = tf(nump, denp);     % modelo contínuo Gp(s) = 1/(10s+1)

Ts = 1;                    % período de amostragem 
Gd  = c2d(Gp, Ts, 'zoh');  % modelo discreto
[numd, dend] = tfdata(Gd, 'v');  
d = 2;                     % atraso 

%% Matriz dinâmica g (resposta ao degrau)
gx(1:5) = 0;           
u(1:5)  = 0;
u(3:N+2) = 1;          
for i = 5:N+3
    gx(i) = -dend(2)*gx(i-1) + numd(2)*u(i-d-1);
end
g = gx(4:length(gx));                     
g = [g  g(length(g))*ones(1, Ny)];       

G = zeros(Ny-N1+1, Nu);
for j = 1:Nu
    for i = 1:Ny-N1+1
        if (N1+i-1-j+1) > 0
            G(i,j) = g(N1+i-1-j+1);
        else
            G(i,j) = 0;
        end
    end
end

%% Malha de simulação
for t = 60:length(yr)-Ny

    % saída da planta
    y(t) = -dend(2)*y(t-1) + numd(2)*u(t-d-1);
    for j = N1:Ny
        cs = zeros(1, N);
        for i = 1:N
            cs(i) = ( g(j+i) - g(i) ) * ( u(t-i) - u(t-i-1) );
        end
        p(j) = y(t) + sum(cs);
        Eo(j-N1+1, 1) = yr(t+j) - p(j);   % erro de predição
    end

    % lei de controle
    DU = (G'*G + lambda*eye(Nu)) \ (G'*Eo);  
    u(t) = u(t-1) + DU(1);                  
end

%% Resultados
t = 1:length(yr)-Ny;
subplot(211), plot(t,y(t), t, yr(t)), ylabel('y(t)'), xlabel('amostra');
subplot(212), plot(t,u(t)), ylabel('u(t)'), xlabel('amostra');
