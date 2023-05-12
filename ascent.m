function dy=ascent(t_0,y_0)
global T g Isp alfa_T h_LP tb S_rocket cD_rocket tT v_w
dy=zeros(5,1);
x=y_0(1);
z=y_0(2);
vx=y_0(3);
vz=y_0(4);
m=y_0(5);
[Te, a, P, rho] = atmosisa(h_LP+z); % thermodynamic conditions at the launch site

Sref=S_rocket;
cD=cD_rocket;

% calculation of angle of vel and apparent vel in the air
if vx==0 && vz==0  % a t=0
    alfa=pi/2;
    V=sqrt(vx^2+vz^2);
    D=0.5*rho*V^2*Sref*cD;
elseif t_0<tb  % during burn vel ~ apparent vel
    alfa=atan2(vz,vx);
    V=sqrt(vx^2+vz^2);
    D=0.5*rho*V^2*Sref*cD;
else  % after burnout
    vx_app=vx-v_w;
    alfa=atan2(vz,vx_app);
    V=sqrt(vx_app^2+vz^2);
    D=0.5*rho*V^2*Sref*cD;
end

%calculation of instantaneous thrust
if t_0>tb
    Tact=0;
    mdot=0;
else
    Tact=interp1(tT,T,t_0);
    mdot=-Tact/g/Isp;
end

dy(1)=vx;
dy(2)=vz;
dy(3)=Tact*sin(alfa_T)/m-D/m*cos(alfa);
dy(4)=Tact*cos(alfa_T)/m-g-D/m*sin(alfa);
dy(5)=mdot;


