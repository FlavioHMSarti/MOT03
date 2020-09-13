clear
close
clc

%% Excel
cp_a = @(Temp) 0.0002*Temp + 0.9449;
cv_a = @(Temp) 0.0002*Temp + 0.6566;
gamma_a = @(Temp) -9/10^5*Temp + 1.4284;
ra = 288.2801664;

cp_f = @(Temp) 0.0003*Temp + 0.936;
cv_f = @(Temp) 0.0003*Temp + 0.6585;
gamma_f = @(Temp) -1/10^4*Temp + 1.4006;
rf = 277.4513479;

cp_p = @(Temp) -6/10^8*Temp^2 + 0.0003*Temp + 0.9508;
cv_p = @(Temp) -6/10^8*Temp^2 + 0.0003*Temp + 0.6638;
gamma_p = @(Temp) 2/10^8*Temp^2 - 0.0001*Temp + 1.4104;
rp = 287.0473442;



%% Hyp
T = zeros(1,8); p = T; V = T;

p_atm = 1;
Tref = 50+273.15;              % on a l'air
COF = 0.8;
Dp_ech = 10/1000;

% F = phi/psi_s
F = 0.7/14.32288;

%% Moteur

D = 79.5/1000;
C = 80.5/1000;
taux_comp = 21.8; % moteur 4b

n_comp = 0.8;
pi_comp = 1.8;
pi_turb = 1.4; 

%% Point 2

T(1) = Tref;
p(1) = p_atm;
p(2) = pi_comp * p(1);
gamma = gamma_a(T(1)); 
T2_is = T(1)*(p(2)/p(1))^((gamma-1)/gamma);

DT = 100;
tol = 1;                %DT = 1K
i = 0;

while (DT > tol && i < 10000)
   Tmoy = (T2_is + T(1))/2;  
   gamma = gamma_a(Tmoy); 
   T2_is = T(1)*(p(2)/p(1))^((gamma-1)/gamma); 
    i = i+1;
    DT = abs(T2_is - 2*Tmoy + T(1));
end
T(2) = T(1) + (T2_is-T(1))/n_comp;

%% Point 3

T(3) = T(2) + COF*(Tref-T(2));
p(3)=p(2)-Dp_ech;

%% Point 3'

% 3' -> dans le moteur, soupape ferm?
r3p = ra;

Cu = pi*D^2*C/4;
V_max = taux_comp * Cu / (taux_comp - 1); % Point mort bas
V_min = Cu/(taux_comp-1);

T3p = T(3);
p3p = p(3); % pas de perte de charge apres le soupape
m3p = p3p*10^5*V_max/r3p/T3p;

%% Point 4
% compression poly 3' - 4
% pv^n = cste

points = csvread('points.csv');
V_cycle = points(:,1);
p_cycle = points(:,2);
%plot(V_cycle,p_cycle)

V_ln = log(V_cycle(3:19));
p_ln = log(p_cycle(3:19));
p1 = polyfit(V_ln,p_ln,1);

n_comp = 1.34; % Excel
p(4) = p3p*taux_comp^n_comp;
T(4) = p(4)*10^5*V_min/m3p/ra;

%% Point 4'
rp = 287.0473442; % Excel

p4p = max(p_cycle);
m4p = m3p;
r4p = rp; % ******
T4p = p4p*10^5 * V_min / m4p / r4p;

%% 4''
p4pp = p4p;
[max_p,indice] = max(p_cycle);
%V4pp = V_cycle(indice);
V4pp = 1.22*V_min;
T4pp = p4pp*10^5 * V4pp / m4p / r4p;

%% 5

p(5)=60; %courbe
V(5)= 3.34*V_min;
T(5) = p(5)*10^5*V(5)/m3p/rp;
%V52 = m3p*rp*T(5)/p(5)/10^5;
DT_45 = abs(T(5)-T4pp)/T4pp*100;

%% 6
n_det = 1.2024; %Excel
p(6) = p(5)*(V(5)/V_max)^n_det;
T(6) = p(6)*10^5*V_max/m3p/rp;
V(6) = m3p*rp*T(6)/p(6)/10^5;

%% 7

DP = abs(6.57964743021319-p(6))/p(6)*100;
p(7) = pi_turb * p_atm;

DT = 100;
tol = 1;                %DT = 1K
i = 0;
gamma = gamma_p(T(6));
T(7) = (inv(gamma)+(gamma-1)/gamma*(p(7)/p(6)))*T(6);

while (DT > tol && i < 10000)
   Tmoy = (T(6) + T(7))/2;
   gamma = gamma_p(Tmoy); 
   T(7) = (inv(gamma)+(gamma-1)/gamma*(p(7)/p(6)))*T(6);
   i = i+1;
   DT = abs(T(7) - 2*Tmoy + T(6));
end


%% 7'
% gaz r?siduels a Padm et Vm
% compresion isentropique des gaz residuels de Pech=P7 a Padm = P3p

p7p = p3p;
DT = 100;
tol = 1;               
i = 0;
gamma = gamma_p(T(7));
T7p = T(7)*(p(7)/p7p)^((1-gamma)/gamma);

while (DT > tol && i < 10000)
   Tmoy = (T(7) + T7p)/2;
   gamma = gamma_p(Tmoy);
   T7p = T(7)*(p(7)/p7p)^((1-gamma)/gamma);
   i = i+1;
   DT = abs(T(7) - 2*Tmoy + T7p);
end

%% Boucle T3p

f = (1/taux_comp)*(T(3)/T(7))*(p(7)/p3p);
tol_DT3p = 1; % 1K
j = 0;

r3p = f*rp+(1-f)*ra; % Diesel
T0 = 273.15;
dif = 100;
T3p_old = T3p;


while (dif > tol_DT3p && j < 10^3)
    j=j+1;
    T03 = (T0+T3p)/2;
    T07 = (T0+T7p)/2;
    cp_03 = f*cp_p(T03)+(1-f)*cp_a(T03);
    T3pn = (f*cp_p(T07)*(T7p-T0)+(1-f)*cp_a(T03)*(T(3)-T0))/cp_03+T0;
    dif = abs(T3pn-T3p_old);
    T3p_old = T3pn;
    m3pn = p3p*10^5*V_max/r3p/T3p;
    f = p(7)*10^5*V_min/m3pn/rp/T(7);
end

%% 3' - 4
Tn = zeros(1,8); pn = T; Vn = T;
r4n = f*rp+(1-f)*ra; %diesel
Tn(4) = p(4)*10^5*V_min/m3pn/r4n;
m4pn = m3pn;

% # 1er BOUCLAGE ####################################################
% m=cst = m3p
%% isochore 4-4'
T4pn = p4p*10^5*V_min/m3pn/rp;
Tmoy = 0.5*(T0+T4pn);
Tmoy2 = 0.5*(T0+Tn(4));
cv = f*cv_p(Tmoy2)+(1-f)*cv_a(Tmoy2);
Qv = m4pn*cv_p(Tmoy)*(T4pn-T0)-m3pn*cv*(Tn(4)-T0);  

%% isobare 4'-4''
m4ppn = m3pn;
T4ppn = p4pp*10^5*V4pp/m4ppn/rp;
Tmoy = 0.5*(T4pn+T4ppn);
Qp = m4ppn * cp_p(Tmoy)*(T4ppn-T4pn);

%% isotherme 4''-5
Tn(5) = p(5)*10^5*V(5)/m3pn/rp;
DT_45 = abs(Tn(5)-T4ppn);
Qt = m3pn * rp * log(V(5)/V4pp);

%% Qapp
Qpp = Qv+Qp+Qt;

%fraction de gaz brules
%xbi = Qi/Qapp
xb4p = Qv/Qpp;
xb4pp = (Qv+Qp)/Qpp;


% # 2eme BOUCLAGE #####################################################3
% m != cst (Diesel)

%% 4-4'
m4pn = (1+F*xb4p)*m3pn;
T4pn = p4p*V_min/m4pn/r3pn;
Tmoy = 0.5*(T0+T4pn);
cv = (1-xb4pp)*cv_a(Tmoy)+xb4p*cv_p(Tmoy);
cv2 = f*cv_p(Tmoy)+(1-f)*cv_a(Tmoy);
Qv = m4pn*cv(Tmoy)*(T4pn-T0)-m3pn*cv2*(Tn(4)-T0); 

%% isobare 4'-4''
m4ppn = (1+F*xb4pp)*m3pn;
r4ppn = (1-xb4pp)*rf+xb4pp*rp;

T4ppn = p4pp*10^5*V4pp/m4ppn/r4ppn;
Tmoy = 0.5*(T4pnn+T0);
Tmoy2 = 0.5*(T4pn+T0);
cp = (1-xb4pp)*cp_f(Tmoy)+xb4pp*cp_p(Tmoy);
cp2 = (1-xb4p)*cp_f(Tmoy)+xb4p*cp_p(Tmoy);
Qp = m4ppn * cp*(T4ppn-T0)-m4pn*cp2*(T4pn-T0);

%% isotherme 4''-5
Vn(5)=V(5);
m5n = (1+F)*m3pn;
r5n = rp;
Tn(5) = p(5)*V(5)/m5n/r5n;
Tmoy = 0.5*(T0+Tn(5));
Tmoy2 = 0.5*(T0+T4ppn);
cv = (1-xb4pp)*cv_f(Tmoy2)+xb4pp*cv_p(Tmoy2);
Qt = m5n*cv_p(Tmoy)*(Tn(5)-T0)-m4ppn*cv*(T4ppn-T0)+0.5*(m4ppn*r4ppn+m5n*r5n)*Tn(5)*log(Vn(5)/V4pp);

Qapp = Qv+Qp+Qt;


% Faire 5-6 et 6-7
