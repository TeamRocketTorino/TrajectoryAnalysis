%% Initialization
clear all;clc;close all
format short
global T g Isp alfa_T h_LP tb S cD tT Sp_body Sp_ogive tp v_w
aVec=[];

%% Input parameters

% ambient
g=9.81; % acceleration of gravity [m/s]
Lat_LP_deg=10; % latitude of the launch site
Lon_LP_deg=30; % longitude of the launch site
v_w = 2.0; % horizontal velocity of wind [m/s]
h_LP=21; % altitude of the launch site [m]

% rocket
m_body_empty = 0.044; % mass of body structure [kg]
m_motor_empty = 0.011; % mass of motor structure [kg]
m_body= m_body_empty + m_motor_empty; % mass of dry body [kg]
m_ogive_empty = 0.078; % mass of ogive structure [kg]
m_payload = 0.022; % mass of payload [kg]
m_ogive=m_ogive_empty+m_payload; % mass of ogive [kg]
mp=0.005; % propellant mass [kg]
m=m_body+m_ogive+mp; % mass of the rocket[kg]
S=pi*((0.035/2)^2); % rocket cross section [m^2]
cD=0.45; % drag coefficient rocket
filename = 'Klima_B4.csv'; % Thrust Curve of the motor
thrust_table = csvread(filename);
tT=thrust_table(:,1); % thrust time [sec]
T=thrust_table(:,2); % thrust values [N]
Itot=5; % total impulse [Ns]
alfa_T_deg=0; % thrust/kick angle [deg]

% parachute
Sp_body=0.094; % cross section body parachute [m^2]
Sp_ogive=0.14; % cross section ogive parachute [m^2]
tp=4; % time delay opening parachute [sec]

%% Processing

% inputs
tb=tT(end); % burning time [sec]
alfa_T=alfa_T_deg/57.3; % thrust/kick angle [rad]
Isp=Itot/(mp*g); % Isp of motor
options = odeset('Events','height','RelTol',1e-13,'AbsTol',1e-15);

% integration

% burn before separation (BURN BABY BURN)
tspan=0:0.002:tb+tp;
y0=[0 0.2 0 0 m];
[t_0, y_0]=ode45(@first_burn,tspan,y0,options);

% body after separation (se separato già un successo)
tspan=tb+tp:0.01:250;
y1=y_0(end,:);
y1(5)=m_body;
[t_1,y_1]=ode45(@body_sep,tspan,y1,options);

% ogive after separation (vedi di non schiantarti)
tspan=tb+tp:0.01:250;
y2=y_0(end,:);
y2(5)=m_ogive;
[t_2,y_2]=ode45(@ogive_sep,tspan,y2,options);

%% Postprocessing

x_0=y_0(:,1); % x position coordinate before separation [m]
z_0=y_0(:,2); % z position coordinate before separation [m]
vx_0=y_0(:,3); % x velocity component before separation [m/s]
vz_0=y_0(:,4); % z velocity component before separation [m/s]
m_0=y_0(:,5); % rocket mass before separation [kg]

x_1=y_1(:,1); % x position coordinate of body after separation [m]
z_1=y_1(:,2); % z position coordinate of body after separation [m]
vx_1=y_1(:,3); % x velocity component of body after separation [m/s]
vz_1=y_1(:,4); % z velocity component of body after separation [m/s]
m_1=y_1(:,5); % rocket mass of body after separation [kg]

x_2=y_2(:,1); % x position coordinate of ogive after separation [m]
z_2=y_2(:,2); % z position coordinate of ogive after separation [m]
vx_2=y_2(:,3); % x velocity component of ogive after separation [m/s]
vz_2=y_2(:,4); % z velocity component of ogive after separation [m/s]
m_2=y_2(:,5); % rocket mass of ogive after separation [kg]

decyear([2021 04 07])

% first burn
for index=1:length(t_0)
    Sref=S;
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

% body after separation
for index=1:length(t_1)
    cD=1.75;
    Sref=Sp_body;
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

% ogive after separation
for index=1:length(t_2)
    cD=1.75;
    Sref=Sp_ogive;
    [Te, a, P, rho] = atmosisa(h_LP+z_2(index));
    if vx_2(index)==0 && vz_2(index)==0 
        alfa=pi/2;
        V=sqrt(vx_2(index)^2+vz_2(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    else
        vx_app=vx_2(index)-v_w;
        alfa=atan2(vz_2(index),vx_app);
        V=sqrt(vx_app^2+vz_2(index)^2);
        D=0.5*rho*V^2*Sref*cD;
    end
    D_list_2(index)=D;
    
    ax_2(index,1)=-D/m_2(index)*cos(alfa);
    
    az_2(index,1)=-g-D/m_2(index)*sin(alfa);
end

%% Representation
% GRAFICI FOR EVERYTHING

figure(1)
plot(x_0,z_0,'-k',x_1,z_1,'-b',x_2,z_2,'-r')
xlabel('x [m]')
ylabel('z [m]')
title('Trajectory')
%axis equal
xlim([-1 4])
ylim([-2 18])
axis equal
grid
legend('ascent','body descent','ogive descent');

figure(2)
plot(t_0,vz_0,'-k',t_1,vz_1,'-b',t_2,vz_2,'-r')
xlabel('t [s]')
ylabel('Vz [m/s]')
ylim([-20 18])
title('Vertical velocity')
grid
legend('ascent','body descent','ogive descent');

figure(4)
plot(t_0,az_0,'-k',t_1,az_1,'-b',t_2,az_2,'-r')
xlabel('t [s]')
ylabel('a_z [m/s^2]')
title('Vertical acceleration')
xlim([0 6]);
ylim([-20 100])
grid
legend('ascent','body descent','ogive descent');

figure(5)
plot(t_0,z_0,'-k',t_1,z_1,'-b',t_2,z_2,'-r')
xlabel('t [s]')
ylabel('z [m]')
title('Motion along z')
grid
legend('ascent','body descent','ogive descent');
ylim([-2 18])

figure(6)
plot(t_0,D_list_0,'-k',t_1,D_list_1,'-b',t_2,D_list_2,'-r')
xlabel('t [s]')
ylabel('D [N]')
ylim([-5 15])
title('Drag')
grid
legend('ascent','body descent','ogive descent');

figure(8)
plot(t_0,Tact_0,'-k')
xlabel('t [s]')
ylabel('T [N]')
title('Thrust')
grid
xlim([0 0.6])