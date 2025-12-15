clear;

%% ===================== PARÁMETROS COMPARTIDOS =====================
T_max = 91;              
Nh = 28180;              
Nm0 = Nh * 5;            

a = 0.30;                
beta_hm = 0.55;          
beta_mh = 0.55;          

T_inc_h = 5;
T_inf_h = 5;
T_inc_m = 8;

sigma_h = 1/T_inc_h;
eta_h   = 1/T_inf_h;
sigma_m = 1/T_inc_m;

mu_m = 1/15;

dt = 1;

times = 0:T_max;

%% =============================================================
%   FUNCIÓN DE SIMULACIÓN GENÉRICA
%% =============================================================
simulate_model = @(Sh,Eh,Ih,Rh,Sm,Em,Im) ...
    simulate(Sh,Eh,Ih,Rh,Sm,Em,Im, ...
             T_max,Nh,a,beta_mh,beta_hm,sigma_h,eta_h,sigma_m,mu_m,dt);


%% =============================================================
%   MODELO A: infección inicial en MOSQUITOS
%% =============================================================
init_Im = 1;
init_Ih = 0;

ShA = Nh;
EhA = 0;
IhA = init_Ih;
RhA = 0;

SmA = Nm0 - init_Im;
EmA = 0;
ImA = init_Im;

[ShA_t,EhA_t,IhA_t,RhA_t,SmA_t,EmA_t,ImA_t] = simulate_model(ShA,EhA,IhA,RhA,SmA,EmA,ImA);


%% =============================================================
%   MODELO B: infección inicial en HUMANOS 
%% =============================================================
init_Ih_B = 1;

ShB = Nh - init_Ih_B;
EhB = 0;
IhB = init_Ih_B;
RhB = 0;

SmB = Nm0;
EmB = 0;
ImB = 0;

[ShB_t,EhB_t,IhB_t,RhB_t,SmB_t,EmB_t,ImB_t] = simulate_model(ShB,EhB,IhB,RhB,SmB,EmB,ImB);


%% =============================================================
%   GRÁFICAS COMPARATIVAS
%% =============================================================

figure;
plot(times,IhA_t,'LineWidth',2); hold on;
plot(times,IhB_t,'LineWidth',2);
legend("Modelo A: infección con mosquitos","Modelo B: infección con humanos");
title("Comparación de Humanos Infectados");
xlabel("Días"); ylabel("Individuos");

figure;
plot(times,ImA_t,'LineWidth',2); hold on;
plot(times,ImB_t,'LineWidth',2);
legend("Modelo A: infección con mosquitos","Modelo B: infección con humanos");
title("Comparación de Mosquitos Infectados");
xlabel("Días"); ylabel("Mosquitos");

figure;
subplot(1,2,1);
plot(times,ShA_t, times, EhA_t, times, IhA_t, times, RhA_t,'LineWidth',2);
legend("S","E","I","R");
title("SEIR Humanos - Modelo A");

subplot(1,2,2);
plot(times,ShB_t, times, EhB_t, times, IhB_t, times, RhB_t,'LineWidth',2);
legend("S","E","I","R");
title("SEIR Humanos - Modelo B");


%% =============================================================
%   FUNCIÓN
%% =============================================================
function [Sh_t,Eh_t,Ih_t,Rh_t,Sm_t,Em_t,Im_t] = simulate(...
    Sh,Eh,Ih,Rh,Sm,Em,Im, ...
    T_max,Nh,a,beta_mh,beta_hm,sigma_h,eta_h,sigma_m,mu_m,dt)

    Sh_t = zeros(1,T_max+1);
    Eh_t = zeros(1,T_max+1);
    Ih_t = zeros(1,T_max+1);
    Rh_t = zeros(1,T_max+1);

    Sm_t = zeros(1,T_max+1);
    Em_t = zeros(1,T_max+1);
    Im_t = zeros(1,T_max+1);

    Sh_t(1)=Sh; Eh_t(1)=Eh; Ih_t(1)=Ih; Rh_t(1)=Rh;
    Sm_t(1)=Sm; Em_t(1)=Em; Im_t(1)=Im;

    for t = 1:T_max

        N_m = Sm + Em + Im;

        lambda_h = a * beta_mh * (Im/Nh);
        lambda_m = a * beta_hm * (Ih/Nh);

        dSh = -lambda_h * Sh;
        dEh =  lambda_h * Sh - sigma_h * Eh;
        dIh =  sigma_h * Eh - eta_h * Ih;
        dRh =  eta_h * Ih;

        births = mu_m * N_m;

        dSm = births - lambda_m * Sm - mu_m * Sm;
        dEm = lambda_m * Sm - sigma_m * Em - mu_m * Em;
        dIm = sigma_m * Em - mu_m * Im;

        Sh = Sh + dSh*dt;
        Eh = Eh + dEh*dt;
        Ih = Ih + dIh*dt;
        Rh = Rh + dRh*dt;

        Sm = Sm + dSm*dt;
        Em = Em + dEm*dt;
        Im = Im + dIm*dt;

        Sh_t(t+1)=Sh; Eh_t(t+1)=Eh; Ih_t(t+1)=Ih; Rh_t(t+1)=Rh;
        Sm_t(t+1)=Sm; Em_t(t+1)=Em; Im_t(t+1)=Im;
    end
end
