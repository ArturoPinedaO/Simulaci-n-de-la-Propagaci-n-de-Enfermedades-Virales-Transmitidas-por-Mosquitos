%% PARAMETROS
T_max = 91;                 % días
Nh = 28180;                % población humana
Nm0 = Nh * 5;              % 10 mosquitos por humano

init_Im = 1;              % <<< MOSQUITOS infectados iniciales
init_Ih = 0;                % <<< Humanos NO inician infectados

a = 0.3;                    % tasa de picadura
beta_hm = 0.55;              % humano -> mosquito
beta_mh = 0.55;              % mosquito -> humano

T_inc_h = 5;
T_inf_h = 5;
T_inc_m = 8;

sigma_h = 1/T_inc_h;
eta_h   = 1/T_inf_h;
sigma_m = 1/T_inc_m;

mu_m = 1/15;                % mortalidad mosquito

%% ESTADOS INICIALES (deterministas)
Sh = Nh;                    % <<< Todos humanos susceptibles
Eh = 0;
Ih = 0;
Rh = 0;

Sm = Nm0 - init_Im;         % susceptibles = total - infectados
Em = 0;
Im = init_Im;               % <<< Se agrega infección inicial de mosquitos

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
dt = 1;

for t = 1:T_max

    N_m = Sm + Em + Im;

    %% Fuerzas de infección
    lambda_h = a * beta_mh * (Im/Nh);   % mosquitos -> humanos
    lambda_m = a * beta_hm * (Ih/Nh);   % humanos -> mosquitos (casi 0 al inicio)

    %% HUMANOS SEIR
    dSh = -lambda_h * Sh;
    dEh =  lambda_h * Sh - sigma_h * Eh;
    dIh =  sigma_h * Eh - eta_h * Ih;
    dRh =  eta_h * Ih;

    %% MOSQUITOS SEI
    births = mu_m * N_m;

    dSm = births - lambda_m * Sm - mu_m * Sm;
    dEm = lambda_m * Sm - sigma_m * Em - mu_m * Em;
    dIm = sigma_m * Em - mu_m * Im;

    %% ACTUALIZAR
    Sh = Sh + dSh;
    Eh = Eh + dEh;
    Ih = Ih + dIh;
    Rh = Rh + dRh;

    Sm = Sm + dSm;
    Em = Em + dEm;
    Im = Im + dIm;

    %% GUARDAR
    Sh_t(t+1)=Sh; Eh_t(t+1)=Eh; Ih_t(t+1)=Ih; Rh_t(t+1)=Rh;
    Sm_t(t+1)=Sm; Em_t(t+1)=Em; Im_t(t+1)=Im;

end

%% GRAFIICAS
times = 0:T_max;

figure; plot(times,Ih_t,'LineWidth',2);
title("Humanos Infectados"); xlabel("Días"); ylabel("Individuos");

figure; plot(times,Im_t,'LineWidth',2);
title("Mosquitos Infectados"); xlabel("Días"); ylabel("Mosquitos");

figure;
plot(times,Sh_t,'LineWidth',2); hold on;
plot(times,Eh_t,'LineWidth',2);
plot(times,Ih_t,'LineWidth',2);
plot(times,Rh_t,'LineWidth',2);
legend("S","E","I","R");
title("SEIR Humanos"); xlabel("Días");
