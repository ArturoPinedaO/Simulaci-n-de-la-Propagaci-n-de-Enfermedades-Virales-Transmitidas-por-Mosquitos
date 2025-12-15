clearvars; close all; rng(42);
%% PARAMETROS
T_max = 91;                 % días
Nh = 28180;                 % población humana
Nm0 = Nh * 5;               % 10 mosquitos por humano

init_Ih = 1;                 % humanos infectados iniciales

a = 0.30;                    % tasa de picadura
beta_hm = 0.55;              % humano -> mosquito
beta_mh = 0.55;              % mosquito -> humano

T_inc_h = 5;
T_inf_h = 5;
T_inc_m = 8;

sigma_h = 1/T_inc_h;
eta_h   = 1/T_inf_h;
sigma_m = 1/T_inc_m;

mu_m = 1/15;                 % mortalidad mosquito


%% ESTADOS INICIALES (deterministas)
Sh = Nh - init_Ih;
Eh = 0;
Ih = init_Ih;
Rh = 0;

Sm = Nm0;
Em = 0;
Im = 0;

%% REGISTROS
Sh_t = zeros(1,T_max+1);
Eh_t = zeros(1,T_max+1);
Ih_t = zeros(1,T_max+1);
Rh_t = zeros(1,T_max+1);

Sm_t = zeros(1,T_max+1);
Em_t = zeros(1,T_max+1);
Im_t = zeros(1,T_max+1);

% guardar estado inicial
Sh_t(1)=Sh; Eh_t(1)=Eh; Ih_t(1)=Ih; Rh_t(1)=Rh;
Sm_t(1)=Sm; Em_t(1)=Em; Im_t(1)=Im;

%% SIMULACIÓN (Euler)
dt = 1;   % paso 1 día

for t = 1:T_max

    N_m = Sm + Em + Im;

    %% Fuerzas de infección
    lambda_h = a * beta_mh * (Im/Nh);   % mosquitos -> humanos
    lambda_m = a * beta_hm * (Ih/Nh);   % humanos -> mosquitos

    %% HUMANOS SEIR
    dSh = -lambda_h * Sh;
    dEh =  lambda_h * Sh - sigma_h * Eh;
    dIh =  sigma_h * Eh - eta_h * Ih;
    dRh =  eta_h * Ih;

    %% MOSQUITOS SEI
    births = mu_m * N_m;   % reposición para mantener población estable

    dSm = births - lambda_m * Sm - mu_m * Sm;
    dEm = lambda_m * Sm - sigma_m * Em - mu_m * Em;
    dIm = sigma_m * Em - mu_m * Im;

    %% ACTUALIZAR
    Sh = Sh + dSh*dt;
    Eh = Eh + dEh*dt;
    Ih = Ih + dIh*dt;
    Rh = Rh + dRh*dt;

    Sm = Sm + dSm*dt;
    Em = Em + dEm*dt;
    Im = Im + dIm*dt;

    %% GUARDAR
    Sh_t(t+1)=Sh; Eh_t(t+1)=Eh; Ih_t(t+1)=Ih; Rh_t(t+1)=Rh;
    Sm_t(t+1)=Sm; Em_t(t+1)=Em; Im_t(t+1)=Im;

end

%% GRAFIICAS
times = 0:T_max;

figure; plot(times,Ih_t,'LineWidth',2);
title("Humanos Infectados (I_h)"); xlabel("Días"); ylabel("Individuos");

figure; plot(times,Im_t,'LineWidth',2);
title("Mosquitos Infectados (I_m)"); xlabel("Días"); ylabel("Mosquitos");

figure;
plot(times,Sh_t,'LineWidth',2); hold on;
plot(times,Eh_t,'LineWidth',2);
plot(times,Ih_t,'LineWidth',2);
plot(times,Rh_t,'LineWidth',2);
legend("S","E","I","R");
title("SEIR Humanos"); xlabel("Días");