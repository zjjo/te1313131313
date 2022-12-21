function [A]=A_calculation(input,Q_nu) 
format long
% h = animatedline();
format long
flag_error=0;
load p_w_RH0_5.mat
NOE_cd=200;
NOE_ad=200;
%% numerical scheme
NOE_wall=5;
difference=1;

% grids for air and coolant
T_air=zeros(NOE_cd,NOE_ad+1);
T_coolant=zeros(NOE_wall+1,NOE_ad);

%% properties
% Air
RH=0.5;
P_air=101325; %[Pa]
vol_air=input(2); %[m/s]
rho=1.225; %[kg/m^2]
CA=input(6)*input(9)/(NOE_cd); %[m^2];
mdot_air = rho*CA*vol_air;  % [kg/s]
AHTC = input(7) / 1000;  %[KW/m^2-K] 
T_air(:,1)=input(3); % K


C_p_air=py.CoolProp.CoolProp.HAPropsSI("cp", "T", T_air(1,1), "R", RH, "P", P_air)/1000; % [KJ/Kg-K]


% R32
C_p_ref=inf;
%C_p_ref = 4193 /1000;  % [KJ/Kg-K]
mdot_ref = input(1);  % [kg/s]
RHTC = input(8) /1000;  % [KW/m^2-K]
T_coolant(1,:) = input(4); % K

% fin
length_ad_fin = input(5);  % [m]
length_cd_fin = input(6); % [m]
t_fin=0.08*0.001;  % [m]

% wall
k_wall = 202.4 / 1000;  % [KW/m-K]
t_wall = 0.5*0.001;  % [m]
length_ad_wall=input(5);
length_cd_wall=input(9)*2+t_fin;

%guess value
intital_guess= (T_air(1,1)+T_coolant(1,1))/2; % K
T_fin =intital_guess*ones(NOE_cd,NOE_ad);
T_wall=intital_guess*ones(NOE_wall,NOE_ad);

% geometry
delta_ad_fin=length_ad_fin/NOE_ad;
delta_cd_fin=length_cd_fin/NOE_cd;
delta_ad_wall=length_ad_wall/NOE_ad;
delta_cd_wall=length_cd_wall/NOE_wall;
delta_A_fin=delta_ad_fin*delta_cd_fin;
delta_A_wall=delta_ad_wall*delta_cd_wall;


%% validation, NTU method
C_p_air_lum=py.CoolProp.CoolProp.HAPropsSI("cp", "T", T_air(1,1), "R", 0.5, "P", 101325)/1000;

C_min=mdot_air*(2*NOE_cd)*C_p_air_lum;
Delta_T_max=abs(T_air(1,1)-T_coolant(1,1));
A_fin_total=2*length_cd_fin*length_ad_fin;
A_wall_total=4*delta_cd_wall*length_ad_fin;
m=sqrt(2*AHTC/(k_wall*t_fin));
eta_dry_rec=tanh(0.5*m*length_cd_fin)/(0.5*m*length_cd_fin);

r1=1/RHTC;
r2=t_fin/(k_wall);

A_total=(A_fin_total+A_wall_total);

surface_eta= 1 - ((A_fin_total / A_total) * (1 - eta_dry_rec));

r3=1/(surface_eta*AHTC);

b=A_total/A_wall_total;

U=1/((r1+r2)*b+r3);

NTU=U * A_total / (mdot_air * (2*NOE_cd) * C_p_air_lum);
epsilon=1-exp(-NTU);
Q_max=C_min*Delta_T_max;
Q_a=epsilon*Q_max;
T_air_out_NTU=T_air(1,1)-Q_a/(mdot_air * (2*NOE_cd) * C_p_air_lum)-273.15;


%% output
A=log(log(1-(Q_nu/Q_a)*(1-exp(-NTU))).^(-1))/log(NTU);
end
