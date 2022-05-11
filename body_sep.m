function dy=body_sep(t_1,y_1)
global T g Isp alfa_T h_LP tb S cD tT Sp_body Sp_ogive tp v_w
dy=zeros(5,1);
x=y_1(1);
z=y_1(2);
vx=y_1(3);
vz=y_1(4);
m=y_1(5);
[Te, a, P, rho] = atmosisa(h_LP+z); % thermodynamic conditions at the launch site

Sref=Sp_body;  % cross section of the body parashute
cD=0.8;  % drag coefficient of the body parashute

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


