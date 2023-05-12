function dy=descent(t_1,y_1)
global g alfa_T h_LP S_parachute cD v_w
dy=zeros(5,1);
x=y_1(1);
z=y_1(2);
vx=y_1(3);
vz=y_1(4);
m=y_1(5);
[Te, a, P, rho] = atmosisa(h_LP+z); % thermodynamic conditions at the launch site

Sref=S_parachute;  % cross section of parashute
cD=0.8;  % drag coefficient of parashute

if vx==0 && vz==0  % a t=0
    alfa=pi/2;
    V=sqrt(vx^2+vz^2);
    D=0.5*rho*V^2*Sref*cD;
else
    vx_app=vx-v_w;
    alfa=atan2(vz,vx_app);
    V=sqrt(vx_app^2+vz^2);
    D=0.5*rho*V^2*Sref*cD;
end

Tact=0;
mdot=0;

dy(1)=vx;
dy(2)=vz;
dy(3)=Tact*sin(alfa_T)/m-D/m*cos(alfa);
dy(4)=Tact*cos(alfa_T)/m-g-D/m*sin(alfa);
dy(5)=mdot;


