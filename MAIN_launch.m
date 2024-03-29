%% Initialization
clear all;clc;close all
format short
global T g Isp alfa_T h_LP tb S_rocket cD_rocket cD_parachute tT S_parachute tp v_w
aVec=[];

%% Input parameters

% ambient
g=9.81; % acceleration of gravity [m/s]
Lat_LP_deg=10; % latitude of the launch site
Lon_LP_deg=30; % longitude of the launch site
v_w = 3; % horizontal velocity of wind [m/s]
h_LP=21; % altitude of the launch site [m]

% rocket
m_body_empty = 0.044; % mass of body structure [kg]
m_motor_empty = 0.011; % mass of motor structure [kg]
m_body= m_body_empty + m_motor_empty; % mass of dry body [kg]
m_ogive_empty = 0.078+0.0900; % mass of ogive structure [kg]
m_payload = 0.022; % mass of payload [kg]
m_ogive=m_ogive_empty+m_payload; % mass of ogive [kg]
mp=0.016; % propellant mass [kg]
m=m_body+m_ogive+mp; % mass of the rocket[kg]
S_rocket=pi*((0.035/2)^2); % rocket cross section [m^2]
cD_rocket=0.4; % drag coefficient rocket
filename = 'Klima_D9.csv'; % Thrust Curve of the motor
thrust_table = csvread(filename);
tT=thrust_table(:,1); % thrust time [sec]
T=thrust_table(:,2); % thrust values [N]
Itot=20.0; % total impulse [Ns]
alfa_T_deg=5; % thrust/kick angle [deg]

% parachute
S_parachute=0.094; % cross section parachute [m^2]
cD_parachute=1.7; % drag coefficient parachute CHANGE CD!!!
tp=3; % time delay opening parachute [sec]

%% Processing

% inputs
tb=2.1; % burning time [sec]
alfa_T=alfa_T_deg/57.3; % thrust/kick angle [rad]
Isp=Itot/(mp*g); % Isp of motor
options = odeset('Events','height','RelTol',1e-13,'AbsTol',1e-15);

% integration

% ascent before separation (BURN BABY BURN)
tspan=0:0.005:tb+tp;
y0=[0 0.4 0 0 m];
[t_0, y_0]=ode45(@ascent,tspan,y0,options);

% descent after separation (se separato già un successo)
tspan=tb+tp:0.01:250;
y1=y_0(end,:);
[t_1,y_1]=ode45(@descent,tspan,y1,options);

%% Postprocessing

x_0=y_0(:,1); % x position coordinate before separation [m]
z_0=y_0(:,2); % z position coordinate before separation [m]
vx_0=y_0(:,3); % x velocity component before separation [m/s]
vz_0=y_0(:,4); % z velocity component before separation [m/s]
m_0=y_0(:,5); % rocket mass before separation [kg]

x_1=y_1(:,1); % x position coordinate of rocket after separation [m]
z_1=y_1(:,2); % z position coordinate of rocket after separation [m]
vx_1=y_1(:,3); % x velocity component of rocket after separation [m/s]
vz_1=y_1(:,4); % z velocity component of rocket after separation [m/s]
m_1=y_1(:,5); % rocket mass of rocket after separation [kg]

decyear([2023 05 30])

% ascent before separation
for index=1:length(t_0)
    Sref=S_rocket;
    cD=cD_rocket;
    [Te, a, P, rho] = atmosisa(h_LP+z_0(index)); % thermodynamic conditions at the launch site
    if vx_0(index)==0 && vz_0(index)==0 
        alfa=pi/2;
        V=sqrt(vx_0(index)^2+vz_0(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    elseif t_0(index)<tb
        alfa=atan2(vz_0(index),vx_0(index));
        V=sqrt(vx_0(index)^2+vz_0(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    else
        vx_app=vx_0(index)-v_w;
        alfa=atan2(vz_0(index),vx_app);
        V=sqrt(vx_app^2+vz_0(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    end
    D_list_0(index)=D;

    if t_0(index)>tb
        Tact_0(index,1)=0;
    else
        Tact_0(index,1)=interp1(tT,T,t_0(index));
    end
    ax_0(index,1)=Tact_0(index,1)*sin(alfa_T)/m_0(index)-D/m_0(index)*cos(alfa);
    
    az_0(index,1)=Tact_0(index,1)*cos(alfa_T)/m_0(index)-g-D/m_0(index)*sin(alfa);
end

% descent after separation
for index=1:length(t_1)
    cD=cD_parachute; 
    Sref=S_parachute;
    [Te, a, P, rho] = atmosisa(h_LP+z_1(index));
    if vx_1(index)==0 && vz_1(index)==0  
        alfa=pi/2;
        V=sqrt(vx_1(index)^2+vz_1(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    else
        vx_app=vx_1(index)-v_w;
        alfa=atan2(vz_1(index),vx_app);
        V=sqrt(vx_app^2+vz_1(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    end
    D_list_1(index)=D;
    
    ax_1(index,1)=-D/m_1(index)*cos(alfa);
    
    az_1(index,1)=-g-D/m_1(index)*sin(alfa);
end

%% Representation
% GRAFICI FOR EVERYTHING

figure(1)
plot(x_0,z_0,'-k',x_1,z_1,'-r')
xlabel('x [m]')
ylabel('z [m]')
title('Trajectory')
%axis equal
%xlim([-1 7])
%ylim([-2 35])
axis equal
grid
legend('ascent','descent');

figure(2)
plot(t_0,vz_0,'-k',t_1,vz_1,'-r')
xlabel('t [s]')
ylabel('Vz [m/s]')
%ylim([-25 25])
title('Vertical velocity')
grid
legend('ascent','descent');

figure(4)
plot(t_0,az_0,'-k',t_1,az_1,'-r')
xlabel('t [s]')
ylabel('a_z [m/s^2]')
title('Vertical acceleration')
%xlim([0 7]);
%ylim([-20 100])
grid
legend('ascent','descent');

figure(5)
plot(t_0,z_0,'-k',t_1,z_1,'-r')
xlabel('t [s]')
ylabel('z [m]')
title('Motion along z')
grid
legend('ascent','descent');
%ylim([-2 35])

figure(6)
plot(t_0,D_list_0,'-k',t_1,D_list_1,'-r')
xlabel('t [s]')
ylabel('D [N]')
%ylim([-5 15])
title('Drag')
grid
legend('ascent','descent');

figure(8)
plot(t_0,Tact_0,'-k')
xlabel('t [s]')
ylabel('T [N]')
title('Thrust')
grid
%xlim([0 1.6])
