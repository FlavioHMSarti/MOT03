clear
close
clc

%% Proprietes physiques
pourc_C = 87.4; % en %
pourc_H = 12.6; % en %
MM = 148.6;    % en g/mol
MM_C = 12;
MM_H = 1;
x = (pourc_C/100)*MM/MM_C;
y = (pourc_H/100)*MM/MM_H;



%% Excel
cp_a = @(Temp) 0.202*Temp + 944.9;
cv_a = @(Temp) 0.202*Temp + 656.59;
gamma_a = @(Temp) -9/10^5*Temp + 1.4284;
ra = 288.2801664;

cp_f = @(Temp) 0.3255*Temp + 935.96;
cv_f = @(Temp) 0.3255*Temp + 658.51;
gamma_f = @(Temp) -1/10^4*Temp + 1.4006;
rf = 277.4513479;

cp_p = @(Temp) -6/10^5*Temp^2 + 0.3304*Temp + 950.81;
cv_p = @(Temp) -6/10^5*Temp^2 + 0.3304*Temp + 663.76;
gamma_p = @(Temp) 2/10^8*Temp^2 - 0.0001*Temp + 1.4104;
rp = 287.0473442;

PCI = 42722.66063; % kJ/kg
phi = 0.7;
psi_s = 14.32288;




%% Hyp
T = zeros(1,8); p = T; V = T;

p_atm = 1;
Tref = 50+273.15;              % on a l'air
COF = 0.8;
Dp_ech = 10/1000;

% F = phi/psi_s
F = phi/psi_s;

%% Moteur

D = 79.5/1000;
C = 80.5/1000;
taux_comp = 21.8; % moteur 4b

n_compresseur = 0.8;
pi_comp = 1.8;
n_turb = 0.8;
pi_turb = 1.4; 

n_meca = 0.875;
n_comb = 0.98;

Ncyl = 4;
N = 3700; %rpm

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
   i = i+1;
   Tmoy = (T2_is + T(1))/2;  
   gamma = gamma_a(Tmoy); 
   T2_is = T(1)*(p(2)/p(1))^((gamma-1)/gamma); 
   DT = abs(T2_is - 2*Tmoy + T(1));
end
T(2) = T(1) + (T2_is-T(1))/n_compresseur;

%% Point 3

T(3) = T(2) + COF*(Tref-T(2));
p(3)=p(2)-Dp_ech;

%% Point 3'

% 3' -> dans le moteur, soupape ferme
r3p = ra;

Cu = pi*D^2*C/4;
V_max = taux_comp * Cu / (taux_comp - 1); % Point mort bas
V_min = Cu/(taux_comp-1);

% Point 3? : Les gaz r?siduels et les gaz frais sont m?lang?s ? la
% pression P3?. On calcule d'abord 3' sans gaz r?siduel et on le
% modifie apr?s dans une boucle
T3p = T(3);
p3p = p(3); % pas de perte de charge apres le soupape
m3p = p3p*10^5*V_max/r3p/T3p;

%% Point 4
% compression poly 3' - 4
% pv^n = cste

points = csvread('Tables\points.csv');
pBas = csvread('Tables\PV_bas.csv'); V_bas = pBas(:,1); p_bas = pBas(:,2);
pHaut = csvread('Tables\PV_haut.csv'); V_haut = pHaut(:,1); p_haut = pHaut(:,2);

n_comp = 1.34; % Excel
p(4) = p3p*taux_comp^n_comp;
T(4) = p(4)*10^5*V_min/m3p/ra;
V(4) = V_min;

%% Point 4'
p4p = max(p_haut);
m4p = m3p;
r4p = rp; % ******
T4p = p4p*10^5 * V_min / m4p / r4p;

%% 4''
p4pp = p4p;
%[max_p,indice] = max(p_cycle);
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
V(7) = (p(6)/p(7))^(1/gamma)*V(6);

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
tol = 1; % 1K
j = 0;

r3p = f*rp+(1-f)*ra; % Diesel
T0 = 273.15;
DT = 100;
T3p_old = T3p;


while (DT > tol && j < 10^3)
    j=j+1;
    T03 = (T0+T3p)/2;
    T07 = (T0+T7p)/2;
    cp_03 = f*cp_p(T03)+(1-f)*cp_a(T03);
    T3pn = (f*cp_p(T07)*(T7p-T0)+(1-f)*cp_a(T03)*(T(3)-T0))/cp_03+T0;
    DT = abs(T3pn-T3p_old);
    T3p_old = T3pn;
    m3pn = p3p*10^5*V_max/r3p/T3p;
    f = p(7)*10^5*V_min/m3pn/rp/T(7);
end


%% 3' - 4
Tn = zeros(1,8); pn = Tn; Vn = Tn;
r4 = f*rp+(1-f)*ra; %diesel
Tn(4) = p(4)*10^5*V_min/m3pn/r4;
m4pn = m3pn;

% # 1er BOUCLAGE ####################################################
% m=cst = m3p
%% isochore 4-4'
T4pn = p4p*10^5*V_min/m4pn/rp;
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
Qt = m3pn * rp * Tn(5) * log(V(5)/V4pp);

%% Qapp
Qapp = Qv+Qp+Qt;

%fraction de gaz brules
%xbi = Qi/Qapp
xv = Qv/Qapp;
xp = (Qv+Qp)/Qapp;


% # 2eme BOUCLAGE #####################################################
% m != cst (Diesel)

%% 4-4'
m4pn = (1+F*xv)*m3pn;
r4p = (1-xv)*rf+xv*rp;
T4pn = p4p*10^5*V_min/m4pn/r4p;
Tmoy = 0.5*(T0+T4pn);
Tmoy2 = 0.5*(T0+Tn(4));
cv = (1-xv)*cv_f(Tmoy)+xv*cv_p(Tmoy);
cv2 = f*cv_p(Tmoy2)+(1-f)*cv_a(Tmoy2);
Qv = m4pn * cv *(T4pn-T0)-m3pn*cv2*(Tn(4)-T0); 

%% isobare 4'-4''
m4ppn = (1+F*xp)*m3pn;
r4pp = (1-xp)*rf+xp*rp;
T4ppn = p4pp*10^5*V4pp/m4ppn/r4pp;

Tmoy = 0.5*(T4ppn+T0);
Tmoy2 = 0.5*(T4pn+T0);
cp = (1-xp)*cp_f(Tmoy)+xp*cp_p(Tmoy);
cp2 = (1-xv)*cp_f(Tmoy2)+xv*cp_p(Tmoy);
Qp = m4ppn * cp * (T4ppn-T0) - m4pn * cp2 * (T4pn-T0);

%% isotherme 4''-5
Vn(5)=V(5);
m5n = (1+F)*m3pn;
r5 = rp;
Tn(5) = p(5)*10^5*V(5)/m5n/r5;

Tmoy = 0.5*(T0+Tn(5));
Tmoy2 = 0.5*(T0+T4ppn);
cv = (1-xp)*cv_f(Tmoy2)+xp*cv_p(Tmoy2);
Qt = m5n*cv_p(Tmoy)*(Tn(5)-T0)-m4ppn*cv*(T4ppn-T0)+0.5*(m4ppn*r4pp+m5n*r5)*Tn(5)*log(Vn(5)/V4pp);
Qappn = Qv+Qp+Qt;

%% 5-6
Tn(6) = p(6)*10^5*V(6)/m5n/rp;
%faire la transformation ou pas

%% 6 - 7

DT = 100;
tol = 1;                %DT = 1K
i = 0;
gamma = gamma_p(Tn(6));
Tn(7) = (inv(gamma)+(gamma-1)/gamma*(p(7)/p(6)))*Tn(6);
Vn(7) = (p(6)/p(7))^(1/gamma)*Vn(6);

while (DT > tol && i < 10000)
   Tmoy = (Tn(6) + Tn(7))/2;
   gamma = gamma_p(Tmoy); 
   Tn(7) = (inv(gamma)+(gamma-1)/gamma*(p(7)/p(6)))*T(6);
   i = i+1;
   DT = abs(Tn(7) - 2*Tmoy + Tn(6));
end


% %% 7' (2) pas besoin?
% 
% DT = 100;
% tol = 1;               
% i = 0;
% gamma = gamma_p(Tn(7));
% T7pn = Tn(7)*(p(7)/p7p)^((1-gamma)/gamma);
% 
% while (DT > tol && i < 10000)
%    Tmoy = (Tn(7) + T7pn)/2;
%    gamma = gamma_p(Tmoy);
%    T7pn = Tn(7)*(p(7)/p7p)^((1-gamma)/gamma);
%    i = i+1;
%    DT = abs(Tn(7) - 2*Tmoy + T7pn);
% end
% %% Boucle T3pn pas besoin?
% 
% %on a f deja
% tol = 1; % 1K
% j = 0;
% 
% r3p = f*rp+(1-f)*ra; % Diesel
% T0 = 273.15;
% DT = 100;
% T3p_old = T3pn;
% fn = f;
% 
% 
% while (DT > tol && j < 10^3)
%     j=j+1;
%     T03 = (T0+T3pn)/2;
%     T07 = (T0+T7pn)/2;
%     cp_03 = fn*cp_p(T03)+(1-fn)*cp_a(T03);
%     T3pn = (fn*cp_p(T07)*(T7pn-T0)+(1-fn)*cp_a(T03)*(T(3)-T0))/cp_03+T0;
%     DT = abs(T3pn-T3p_old);
%     T3p_old = T3pn;
%     m3pn = p3p*10^5*V_max/r3p/T3pn;
%     fn = p(7)*10^5*V_min/m3pn/rp/Tn(7);
% end
% Df = abs(fn-f)/f*100;

Tn(Tn==0)=T(Tn==0);
pn(pn==0)=p(pn==0);
Vn(Vn==0)=V(Vn==0);
% Tn=Tn'; pn=pn'; Vn=Vn';
% T=T'; p=p'; V=V';

%% 7 - 8
tol_q = 15/100; %qwaste/qturb
j = 0;
qwaste = 1; qturb = 1;
pi_turbn = pi_turb;
rapport_q_old = 0;

while(qwaste/qturb > tol_q  && j < 10^5)
    
    if (rapport_q_old > qwaste/qturb)
        pi_turbn = pi_turb*(1-0.1);
    else
        pi_turbn = pi_turb*(1+0.1);
    end

p(8) = p(7)/pi_turbn;

DT = 100;
tol = 1;                %DT = 1K
i = 0;

gamma = gamma_p(Tn(7));
T8_is = Tn(7)*(p(8)/p(7))^((gamma-1)/gamma); 

while (DT > tol && i < 10000)
   i = i+1;
   Tmoy = (T8_is + Tn(7))/2;  
   gamma = gamma_p(Tmoy); 
   T8_is = Tn(7)*(p(8)/p(7))^((gamma-1)/gamma); 
   DT = abs(T8_is - 2*Tmoy + Tn(7));
end
Tn(8) = Tn(7) + (T8_is-Tn(7))/n_turb;

%% Calcul des debits massiques (kg/s)
Tmoy = 0.5*(T(1)+T(2));
Tmoy2 = 0.5*(Tn(7)+Tn(8));

m3 = p3p*10^5*V_max/r3p/T3pn*(1-f);
qcomp = m3 * N/120*Ncyl;
qturb = -qcomp*cp_a(Tmoy)*(T(2)-T(1))/(n_meca*(cp_p(Tmoy2)*(Tn(8)-T(7))));

qcarb = m3p*(1-f)*F*N/120*Ncyl;
rapport_q_old = qwaste/qturb;
qwaste = qcomp+qcarb-qturb;


j = j+1;

end

Dq = (qwaste)/(qcomp+qcarb+qturb+qwaste)*100;

%% 7 - 8'
p8p = p_atm;
T8p = Tn(7);

%% 8' - 9
qmoteur = qwaste + qturb;
Tn(9) = qwaste/qmoteur*T8p + qturb/qmoteur*Tn(8);

%% Calcul des travaux
r3p = (1-f)*ra+f*rp;
W3p4 = m3pn*r3p/(n_comp-1)*(Tn(4)-T3pn);
Wp = -max(p_haut)*10^5*(V4pp-V_min);

mr = 0.5*(m4ppn*r4pp+m5n*rp);
Wt = -mr*Tn(5)*log(Vn(5)/V4pp);

W56 = m5n*rp/(n_det-1)*(Tn(6)-Tn(5));

W_hp = W3p4+Wp+Wt+W56;
W_bp = (pn(7)-p3p)*10^5*(V_max-V_min);
W_th = W_hp + W_bp;

%% Bilan thermique
ma = m3p*(1-f);
Qcalo = ma * F * PCI*1000;
Qcomb = n_comb * Qcalo;
Qapp = Qv+Qp+Qt;

%% Calcul des chaleurs
Tmoy = 0.5*(Tn(2)+Tn(3));
m_air = qcomp / (N/120 * Ncyl);
Q23 = m_air * cp_a(Tmoy) * (Tn(3)-Tn(2)); % < 0

Tmoy = 0.5*(T3pn+Tn(4));
cv = f*cv_p(Tmoy)+(1-f)*cv_a(Tmoy);
gamma = f*gamma_p(Tmoy)+(1-f)*gamma_a(Tmoy);
Q3p4 = m3pn*cv*(Tn(4)-T3pn)-W3p4; % si n_comp = 1.34 < gamma = 1.37, Q < 0

Tmoy = 0.5*(Tn(5)+Tn(6));
cv = f*cv_p(Tmoy)+(1-f)*cv_a(Tmoy);
gamma = f*gamma_p(Tmoy)+(1-f)*gamma_a(Tmoy);
Q56 = m3pn*cv*(Tn(6)-Tn(5))-W56; % si n_det < gamma, Q > 0

Q_pertes_comb = Qcomb-Qcalo; % < 0
Q_pertes_paroi = Qapp - Qcomb; % <0

Tmoy = 0.5*(Tn(7)+Tn(8));
Wturb = qturb/(N/120*Ncyl) * cp_p(Tmoy) * (Tn(8)-Tn(7));
Q_frottements = Wturb*(1-n_meca); % <0

Tmoy =  0.5*(Tn(1)+T0);
Tmoy2 = 0.5*(Tn(9)+T0);
m_echap = qmoteur/(N/120*Ncyl);
Q_echap = m_air * cp_a(Tmoy) * (Tn(1)-T0) - m_echap * cp_p(Tmoy2) * (Tn(9)-T0); % <0

Pertes_moteur = Q_pertes_comb + Q_pertes_paroi + Q3p4 + Q56; % sans echap
Pertes_turbo = Q23 + Q_frottements;

indice = 1; %  1 -> total   2-> Pertes moteur   3-> Pertes turbo

if (indice == 1)
    X = [Pertes_moteur, Pertes_turbo, W_th, Q_echap];
    labels = {'Pertes Moteur', 'Pertes Turbo', 'Travail Recup', 'Q echap'};
elseif (indice ==2)
    X = [Q_pertes_comb, Q_pertes_paroi, Q3p4, Q56, Q_echap];
    labels = {'Q pertes comb', 'Q pertes paroi', 'Q3p4', 'Q56' ,'Q echap'};
else 
    X = [Q23, Q_frottements];
    labels = {'Q23', 'Q frottements'};
end

% pie(abs(X))
% legend(labels,'Location','southeast')


% discuter comment diminuer les pertes

DQ = (Qcalo-abs(Pertes_moteur+Q_echap+Pertes_turbo+W_th))/Qcalo*100;

%% Performances

Pi = abs(W_th) * N/120 * 4; %on a 4 cylindres
PMI = abs(W_th)/(V_max-V_min); %on dit que W_th = Wi, c'est le cycle enveloppee
PMI_hp = abs(W_hp)/(V_max-V_min);
PMI_bp = abs(W_bp)/(V_max-V_min);
CSI = qcarb/Pi;
CSI_g = CSI*(1000*3600*1000);
n_icycle = abs(W_th)/Qcalo;

% 
% n_th = -W_th / Qcomb;
% n_cycle = W_i / W_th;
% n_m = W_e / W_i; %mecanique de l'ensemble
% n_i = n_comb * n_th * n_cycle;
% n_eff = n_i * n_m;


%% Temperature adiabatique de flamme 
% prevoir les oxydes
% vitesse de formation des oxydes

phi1=1; % le melange air/carburant est a richesse 1 (meme pour diesel)

%% a) A partir du PCI
qp = phi1*PCI*1000/(phi1+psi_s);
T5a = Tn(4) + qp/cp_p(Tn(4));

DT = 100;
tol = 1;                %DT = 1K
i = 0;

while (DT > tol && i < 10000)
   i = i+1;
   Tmoy = 0.5*(Tn(4)+T5a);
   T5a = Tn(4) + qp/cp_p(Tmoy);
   DT = abs(T5a-2*Tmoy+Tn(4));
end


%% A partir de la chaleur de reaction
tCxHy = csvread('Tables\CxHy_Tables.csv');
tO2 = csvread('Tables\O2_Tables.csv');
tN2 = csvread('Tables\N2_Tables.csv');
tH2O = csvread('Tables\H2O_Tables.csv');
tCO2 = csvread('Tables\CO2_Tables.csv');

% combustion
    % reactifs
H_C  = 1*interp_data(tCxHy,Tn(4),1,2);
H_O2 = (x+y/4)*interp_data(tO2,Tn(4),1,3);
H_N2 = (x+y/4)*3.76*interp_data(tN2,Tn(4),1,3);
    % produits
H_CO2  = x*interp_data(tCO2,Tn(4),1,3);
H_H2O  = y/2*interp_data(tH2O,Tn(4),1,3);
H_N2_2 = (x+y/4)*3.76*interp_data(tN2,Tn(4),1,3);
H_produits = (H_CO2+H_H2O+H_N2_2);
H_reactifs = (H_C+H_O2+H_N2);
H1 = H_produits - H_reactifs; H1 = abs(H1);

DH = 1000;
tol = 1/100;
i = 0;
T5b = Tn(4)+50;

while (DH > tol && i < 10^4)
   i = i+1;
   % echauffement % Hs(T5)-Hs(T4)
   H_CO2 = x*(interp_data(tCO2,T5b,1,2)-interp_data(tCO2,Tn(4),1,2));
   H_H2O = y/2*(interp_data(tH2O,T5b,1,2)-interp_data(tH2O,Tn(4),1,2));
   H_N2_2 = (x+y/4)*3.76*(interp_data(tN2,T5b,1,2)-interp_data(tN2,Tn(4),1,2));
   H2 = H_CO2+H_H2O+H_N2_2;
   DH = abs(H1-H2)/H1;
   T5b = T5b+50;
end

DT_ver = abs(T5a-T5b);

%% Plot
% Vplot = [V_max; V_min; V_min; V4pp; Vn(5); V_max; V_max];
% pplot = [p3p; pn(4); p4p; p4pp; pn(5); pn(6); p3p];


% loglog(V_bas*V_min,p_bas,V_haut*V_min,p_haut, Vplot, pplot, '*')
% loglog(Vplot,pplot, '*')


