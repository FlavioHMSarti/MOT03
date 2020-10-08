proprietes_gaz

%% Proprietes gaz

load('input_variables/proprietes_gaz.mat')
F = phi/psi_s;


% Rapport



%% Moteur

T = zeros(1,8); p = T; V = T;

moteur = readtable('input_variables/moteur_spec.csv','HeaderLines',1);
moteur = table2array(moteur(:,2));

D = moteur(1);
C = moteur(2);
taux_comp = moteur(3);
n_compresseur = moteur(4);
pi_comp = moteur(5);
n_turb = moteur(6);
pi_turb = moteur(7); 
n_meca = moteur(8);
n_comb = moteur(9);
Ncyl = moteur(10);
N = moteur(11);

T_capot = moteur(12) + 273.15;
T_ref_ech = moteur(13) + 273.15;
p_atm = moteur(14);
COF = moteur(15);
Dp_ech = moteur(16)/1000; % mbar - bar



%% Point 2

T(1) = T_capot;
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

T(3) = T(2) + COF*(T_ref_ech-T(2));
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

pBas = csvread('input_variables\PV_bas.csv',1); V_bas = pBas(:,1); p_bas = pBas(:,2);
pHaut = csvread('input_variables\PV_haut.csv',1); V_haut = pHaut(:,1); p_haut = pHaut(:,2);

% [~,i] = max(V_bas);
% p3p = p_bas(i);

T3p = T(3);
p3p = p(3); % pas de perte de charge apres le soupape
% [~,indice] = max(V_bas);
% p3p = p_bas(indice);
m3p = p3p*10^5*V_max/r3p/T3p;

%% Point 4
% compression poly 3' - 4
% pv^n = cste

n_comp = polytropique(p_bas,V_bas,p3p);
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
V4pp = 1.22*V_min;
T4pp = p4pp*10^5 * V4pp / m4p / r4p;

%% 5

[p(5),V(5),indices] = trouver_isotherme(p_haut, V_haut, p4pp, V4pp/V_min);
V(5) = V(5) * V_min;

T(5) = p(5)*10^5*V(5)/m3p/rp;
%V52 = m3p*rp*T(5)/p(5)/10^5;
DT_45 = abs(T(5)-T4pp)/T4pp*100;

%% 6
n_det = polytropique(p_haut,V_haut,p(5));
p(6) = p(5)*(V(5)/V_max)^n_det;
T(6) = p(6)*10^5*V_max/m3p/rp;
V(6) = m3p*rp*T(6)/p(6)/10^5;

%% 7

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

Tn(Tn==0)=T(Tn==0);
pn(pn==0)=p(pn==0);
Vn(Vn==0)=V(Vn==0);
% *****
p7p = p3p;
DT = 100;
tol = 1;               
i = 0;
gamma = gamma_p(Tn(7));
T7pn = Tn(7)*(pn(7)/p7p)^((1-gamma)/gamma);

while (DT > tol && i < 10000)
   Tmoy = (Tn(7) + T7pn)/2;
   gamma = gamma_p(Tmoy);
   T7pn = T(7)*(pn(7)/p7p)^((1-gamma)/gamma);
   i = i+1;
   DT = abs(Tn(7) - 2*Tmoy + T7pn);
end
fn = (1/taux_comp)*(Tn(3)/Tn(7))*(p(7)/p3p);
tol = 1; % 1K

r3p = f*rp+(1-f)*ra; % Diesel
T0 = 273.15;
T3p_old = T3pn;

DT = 100;
j = 0;
while (DT > tol && j < 10^3)
    j=j+1;
    T03 = (T0+T3pn)/2;
    T07 = (T0+T7p)/2;
    cp_03 = f*cp_p(T03)+(1-fn)*cp_a(T03);
    T3pn = (f*cp_p(T07)*(T7pn-T0)+(1-fn)*cp_a(T03)*(Tn(3)-T0))/cp_03+T0;
    DT = abs(T3pn-T3p_old);
    T3p_old = T3pn;
    m3pn = p3p*10^5*V_max/r3p/T3pn;
    fn = p(7)*10^5*V_min/m3pn/rp/Tn(7);
end

%% 7 - 8
tol_q = 15/100; %qwaste/qturb
j = 0;
qwaste = 100; qturb = 50;
pi_turbn = pi_turb;
rapport_q_old = 0;

while(abs(qwaste)/qturb > tol_q  && j < 10^5)
 
    if (rapport_q_old > qwaste/qturb)
        pi_turbn = pi_turb*(1-0.1);
    else
        pi_turbn = pi_turb*(1+0.1);
    end

pn(8) = pn(7)/pi_turbn;

DT = 100;
tol = 1;                %DT = 1K
i = 0;

gamma = gamma_p(Tn(7));
T8_is = Tn(7)*(pn(8)/p(7))^((gamma-1)/gamma);


while (DT > tol && i < 10000)
   i = i+1;
   Tmoy = (T8_is + Tn(7))/2;  
   gamma = gamma_p(Tmoy); 
   T8_is = Tn(7)*(pn(8)/p(7))^((gamma-1)/gamma); 
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
Dq = abs(qwaste)/(qcomp+qcarb+qturb+abs(qwaste))*100;

%% 7 - 8'
p8p = p_atm;
T8p = Tn(7);

%% 8' - 9
qmoteur = qwaste + qturb;
Tn(9) = qwaste/qmoteur*T8p + qturb/qmoteur*Tn(8);
pn(9) = p_atm;

%% Bilan thermique
ma = m3p*(1-f);
Qcalo = ma * F * PCI*1000;
Qcomb = n_comb * Qcalo;
Qapp = Qv+Qp+Qt;
xvn = Qv/Qapp;
xpn = (Qv+Qp)/Qapp;


%% Calcul des travaux
r3p = (1-f)*ra+f*rp;
W3p4 = m3pn*r3p/(n_comp-1)*(Tn(4)-T3pn);
W4p4pp = -max(p_haut)*10^5*(V4pp-V_min);

mr = 0.5*(m4ppn*r4pp+m5n*rp);
W4pp5 = -mr*Tn(5)*log(Vn(5)/V4pp);

W56 = m5n*rp/(n_det-1)*(Tn(6)-Tn(5));

W_hp = W3p4+W4p4pp+W4pp5+W56;
W_bp = (pn(7)-p3p)*10^5*(V_max-V_min);
W_th = W_hp + W_bp;


Tmoy = 0.5*(T0+T4pn);
Tmoy2 = 0.5*(T0+Tn(4));
cv = (1-xvn)*cv_f(Tmoy)+xvn*cv_p(Tmoy);
cv2 = f*cv_p(Tmoy2)+(1-f)*cv_a(Tmoy2);
Qv = m4pn * cv *(T4pn-T0)-m3pn*cv2*(Tn(4)-T0);

Tmoy = 0.5*(T4ppn+T0);
Tmoy2 = 0.5*(T4pn+T0);
cp = (1-xpn)*cp_f(Tmoy)+xpn*cp_p(Tmoy);
cp2 = (1-xvn)*cp_f(Tmoy2)+xvn*cp_p(Tmoy);
Qp = m4ppn * cp * (T4ppn-T0) - m4pn * cp2 * (T4pn-T0);


Tmoy = 0.5*(T0+Tn(5));
Tmoy2 = 0.5*(T0+T4ppn);
cv = (1-xpn)*cv_f(Tmoy2)+xpn*cv_p(Tmoy2);
Qt = m5n*cv_p(Tmoy)*(Tn(5)-T0)-m4ppn*cv*(T4ppn-T0)+0.5*(m4ppn*r4pp+m5n*r5)*Tn(5)*log(Vn(5)/V4pp);



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
Q_pertes_paroi = Qappn - Qcomb; % <0

Tmoy = 0.5*(Tn(7)+Tn(8));
Wturb = qturb/(N/120*Ncyl) * cp_p(Tmoy) * (Tn(8)-Tn(7));
Q_frottements = Wturb*(1-n_meca); % <0

Tmoy =  0.5*(Tn(1)+T0);
Tmoy2 = 0.5*(Tn(9)+T0);
m_echap = qmoteur/(N/120*Ncyl);
Q_echap = m_air * cp_a(Tmoy) * (Tn(1)-T0) - m_echap * cp_p(Tmoy2) * (Tn(9)-T0); % <0

indice = 2;

if (indice == 1)
    X = [Q_pertes_comb,Q_pertes_paroi,Q_frottements,Q_echap,Q23];
    labels = {'Pertes Combustion', 'Pertes Paroi', 'Frottement', 'Q echap','Q23'};
    
elseif (indice ==2)
    X = [abs(Q_pertes_comb+Q_frottements+Q23), abs(Q_pertes_paroi), abs(Q_echap), abs(W_th)];
    labels = {'Autres pertes', 'Pertes Paroi', 'Q echap', 'Travail Utile Cycle'};

else 
    X = [Q23, Q_frottements];
    labels = {'Q23', 'Q frottements'};
end

if (opt_plot_1)
    figure()
    pie(abs(X))
    legend(labels,'Location','bestoutside','Orientation','vertical')
    title('Distribution energetique du systeme')

end

%% Performances

Pi = abs(W_th) * N/120 * 4; %on a 4 cylindres
PMI = abs(W_th)/(V_max-V_min); %on dit que W_th = Wi, c'est le cycle enveloppee
PMI_hp = abs(W_hp)/(V_max-V_min);
PMI_bp = abs(W_bp)/(V_max-V_min);
CSI = qcarb/Pi;
CSI_g = CSI*(1000*3600*1000);
n_icycle = abs(W_th)/Qcalo;


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
% Les colonnes
% T	 Cp	 s_l	h_l 	h_s	 Dh_f	 DG_f	 logKp

col  = 6;   % Dh_f
col2 = 5;   % h_s

% combustion
    % reactifs
H_C  = 1*interp_data(tCxHy,Tn(4),1,2);
H_O2 = (x+y/4)*interp_data(tO2,Tn(4),1,col);
H_N2 = (x+y/4)*3.76*interp_data(tN2,Tn(4),1,col);
    % produits
H_CO2  = x*interp_data(tCO2,Tn(4),1,col);
H_H2O  = y/2*interp_data(tH2O,Tn(4),1,col);
H_N2_2 = (x+y/4)*3.76*interp_data(tN2,Tn(4),1,col);

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
   H_CO2 = x*(interp_data(tCO2,T5b,1,col2)-interp_data(tCO2,Tn(4),1,col2));
   H_H2O = y/2*(interp_data(tH2O,T5b,1,col2)-interp_data(tH2O,Tn(4),1,col2));
   H_N2_2 = (x+y/4)*3.76*(interp_data(tN2,T5b,1,col2)-interp_data(tN2,Tn(4),1,col2));
   H2 = H_CO2+H_H2O+H_N2_2;
   DH = abs(H1-H2)/H1;
   T5b = T5b+50;
end

DT_ver = abs(T5a-T5b);


%% Plot
if(opt_plot_2)
    figure()
    Vplot = [V_max/V_min;   1;      1;      V4pp/V_min;    Vn(5)/V_min;    V_max/V_min;     V_max/V_min];
    pplot = [p3p;           pn(4);  p4p;    p4pp;          pn(5);          pn(6);           p3p];
    plot(V_bas,p_bas,'blue',V_haut,p_haut,'blue',Vplot, pplot,'black*','LineWidth',2)
    title('Cycle Thermodynamique (Reel vs Enveloppe)','fontweight','bold')
    xlabel('V/V_0','fontweight','bold')
    ylabel('Pression (bar)','fontweight','bold')
    legend('Enveloppe')
    grid()
end

if(opt_plot_3)
    figure()
    Vplot = [V_max/V_min;   1;      1;      V4pp/V_min;    Vn(5)/V_min;    V_max/V_min;     V_max/V_min];
    pplot = [p3p;           pn(4);  p4p;    p4pp;          pn(5);          pn(6);           p3p];
    loglog(Vplot, pplot,'*-black',V_bas,p_bas,'blue',V_haut,p_haut,'blue','LineWidth',2)
    title('Cycle Thermodynamique (Reel vs Enveloppe)','fontweight','bold')
    xlabel('V/V_0','fontweight','bold')
    ylabel('Pression (bar)','fontweight','bold')
    legend('Enveloppe')
    grid()
end


%% Rapport des resultats

outVar = [PCI/1000; n_comp; n_det;f*100 ;fn*100;pi_turbn; qcomp*1000; qturb*1000; ...
    abs(qwaste*1000);qcarb*1000; W3p4*4;W4p4pp*4 ;W4pp5*4;W56*4; W_hp*4;...
    W_bp*4; W_th*4;Wturb*4;Q23*4;Q3p4*4;Qv*4;Qp*4;Qt*4;Q56*4;Q_frottements*4;...
    Qcalo*4; Qcomb*4; Qapp*4; Q_pertes_comb*4;Q_pertes_paroi*4; Q_echap*4;...
    PMI_hp/10^5; PMI_bp/10^5;...
PMI/10^5;Pi/10^3;CSI_g;n_icycle*100;T5a;T5b];
namesVar = ["PCI";"n_comp";"n_det";"f";"fn";"pi_det";"q_comp";"q_turb";"q_waste";...
    "q_carb";"W_3p4'"; "W_4'4''";"W_4''5";"W_56"; "W_hp"; "W_bp";...
    "W_cycle";"W_turbine";"Q_23";"Q_3'4";"Q_44'";"Q_4'4''";"Q_4''5";"Q_56";...
    "Q_frottement";"Q_calo";"Q_comb" ; "Q_app"; "Q_pertes_combustion"; ...
    "Q_pertes_paroi";"Q_echap";...
    "PMI_hp";"PMI_bp";"PMI";"Pi";"CSI"; "n_cycle";"T_5_PCI";"T_5_energie"];
uniVar = ["kJ/kg";"-";"-";"%";"%";"-";"g/s";"g/s";"g/s";...
    "g/s";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"J";"bar";...
    "bar";"bar";"kW";"g/kW-h";"%";"K";"K"];

T1 = table(namesVar,uniVar,outVar);

T2_T = [Tn(1:3)';T3p;Tn(4);T4pn;T4ppn;Tn(5);Tn(6:7)';T7pn;Tn(8);T8p;Tn(9)];
T2_p = [pn(1:3)';p3p;pn(4);p4p;p4pp;pn(5);pn(6:7)';p7p;pn(8);p8p;pn(9)];
T2_names = [ "1 "; "2 "; "3"; "3' "; "4 "; "4' " ; "4'' "; "5 "; "6 "; "7";"7'";"8";"8'";"9"];
T2 = table(T2_names,T2_T,T2_p);

T1.Properties.VariableNames = {'Variable', 'Unite', 'Value'};
T2.Properties.VariableNames = {'Point','T_K', 'p_bar'};

output_folder = 'output_variables';
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end
write(T1,'output_variables/rapport.csv','Delimiter',',')
write(T2,'output_variables/rapport_pT.csv','Delimiter',',')