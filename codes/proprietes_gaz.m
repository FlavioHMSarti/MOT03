% constantes
MM_C = 12/1000;                  % kg/mol       
MM_H = 1/1000;                   % kg/mol
MM_O = 16/1000;                  % kg/mol
MM_N = 14/1000;                  % kg/mol
MM_CO2 = MM_C+2*MM_O;            % kg/mol
MM_H2O = 2*MM_H + MM_O;          % kg/mol
R = 8.314;                       % J/mol-K

% carburant

carburant = readtable('input_variables/carburant_spec.csv','HeaderLines',1);
carburant = table2array(carburant(:,2));

pMasse_C = carburant(1)/100;         % "%"
pMasse_H = carburant(2)/100;         % "%"
pMasse_O = carburant(3)/100;         % "%"
MM_carb  = carburant(4)/1000;        % kg/mol

% moteur
phi = carburant(5);                  % richesse

% reaction
% C_x H_y O_z + (x + y/4 - z/2)*(O2 + 3.78*N2) = 
% x CO2 + y/2 * H2O + 3.78*(x + y/4 - z/2)*N2

x = pMasse_C*MM_carb/MM_C;
y = pMasse_H*MM_carb/MM_H;
z = pMasse_O*MM_carb/MM_O;

% tables thermodynamiques
tCxHy = csvread('Tables\CxHy_Tables.csv');
tO2 = csvread('Tables\O2_Tables.csv');
% Les colonnes
% T	 Cp	 s_l	h_l 	h_s	 Dh_f	 DG_f	 logKp

tN2 = csvread('Tables\N2_Tables.csv');
tH2O = csvread('Tables\H2O_Tables.csv');
tCO2 = csvread('Tables\CO2_Tables.csv');

%% air
pVol_O = 21/100;                                    % "%"
pVol_N = 79/100;                                    % "%"
MM_air = (pVol_O*(MM_O*2)+pVol_N*(MM_N*2));         % kg/mol
ra = R / MM_air;                                	% J/kg-K   

T = [298.15, 500, 1000, 1500, 2000, 2500, 3000]; cp_air = zeros(1,3);
for i=1:length(T)
    cp_air(i) = interp_data(tN2,T(i),1,2)*pVol_N;
    cp_air(i) = cp_air(i) + interp_data(tO2,T(i),1,2)*pVol_O;
end

cp_air = cp_air/MM_air;                             % J/kg-K
cv_air = cp_air - ra;                               % J/kg-K
gamma_air = cp_air./(cp_air-ra);

%% gaz frais
MM_frais = (phi*MM_carb + (x+y/4)*2*MM_O + 3.76*(x+y/4)*2*MM_N)/(phi+(x+y/4)*(1+3.76));     % kg/mol
rf = R / MM_frais;                                                                     % J/kg-K

cp_frais = zeros(1,3);
for i=1:length(T)
    cp_frais(i) = phi * interp_data(tCxHy, T(i),1,4);
    cp_frais(i) = cp_frais(i) + interp_data(tO2,T(i),1,2) * (x+y/4);
    cp_frais(i) = cp_frais(i) + interp_data(tN2,T(i),1,2) * 3.76 * (x+y/4);
    cp_frais(i) = cp_frais(i) / (phi+(x+y/4)*(1+3.76));
end

cp_frais = cp_frais/MM_frais;                         % J/kg-K
cv_frais = cp_frais - rf;                             % J/kg-K
gamma_frais = cp_frais./(cp_frais-rf);

%% gaz brules (produits)

% composition volumique
denominateur = phi*x + phi*(y/2) + (1-phi)*(x+y/4) + (x+y/4)*3.76;
pVol_CO2 = phi*x/denominateur;
pVol_H2O = phi*(y/2)/denominateur;
pVol_O22 = (1-phi)*(x+y/4)/denominateur;
pVol_N22 = 3.76*(x+y/4)/denominateur;

MM_prod = (phi*x*MM_CO2 + phi*(y/2)*MM_H2O + (1-phi)*(x+y/4)*MM_O*2 + (x+y/4)*3.76*MM_N*2)/denominateur;
rp = R/MM_prod;

cp_prod = zeros(1,3);
for i=1:length(T)
    cp_prod(i) = interp_data(tCO2, T(i),1,2) * phi * x;
    cp_prod(i) = cp_prod(i) + interp_data(tH2O,T(i),1,2) * phi * (y/2);
    cp_prod(i) = cp_prod(i) + interp_data(tO2,T(i),1,2) * (1-phi) * (x+y/4);
    cp_prod(i) = cp_prod(i) + interp_data(tN2,T(i),1,2) * 3.76 * (x+y/4);
    cp_prod(i) = cp_prod(i) / denominateur;
end

cp_prod = cp_prod/MM_prod;                     % J/kg-K
cv_prod = cp_prod - rp;                        % J/kg-K
gamma_prod = cp_prod./(cp_prod-rp);

%% interpolation
cp_a = polyfit(T(1:3),cp_air(1:3),1);
cv_a = polyfit(T(1:3),cv_air(1:3),1);
gamma_a = polyfit(T(1:3),gamma_air(1:3),1);

cp_f = polyfit(T(1:3),cp_frais(1:3),1);
cv_f = polyfit(T(1:3),cv_frais(1:3),1);
gamma_f = polyfit(T(1:3),gamma_frais(1:3),1);

cp_p = polyfit(T,cp_prod,2);
cv_p = polyfit(T,cv_prod,2);
gamma_p = polyfit(T,gamma_prod,2);

cp_a    = @(Temp) polyval(cp_a, Temp);
cv_a    = @(Temp) polyval(cv_a, Temp);
gamma_a = @(Temp) polyval(gamma_a, Temp);
cp_f    = @(Temp) polyval(cp_f, Temp);
cv_f    = @(Temp) polyval(cv_f, Temp);
gamma_f = @(Temp) polyval(gamma_f, Temp);
cp_p    = @(Temp) polyval(cp_p, Temp);
cv_p    = @(Temp) polyval(cv_p, Temp);
gamma_p = @(Temp) polyval(gamma_p, Temp);

%% PCI
% Les colonnes
% T	 Cp	 s_l	h_l 	h_s	 Dh_f	 DG_f	 logKp
Tref = 298.15;
Hf_CO2 = interp_data(tCO2, Tref,1,6);
Hf_H2O = interp_data(tH2O, Tref,1,6);
Hf_CxHy = interp_data(tCxHy, Tref,1,2);
Hf_O2 = interp_data(tO2, Tref,1,6);
Hf_N2 = interp_data(tN2, Tref,1,6);

% reaction
% C_x H_y O_z + (x + y/4 - z/2)*(O2 + 3.78*N2) = 
% x CO2 + y/2 * H2O + 3.78*(x + y/4 - z/2)*N2

DH_comb = x*Hf_CO2 + y/2*Hf_H2O + 3.78*(x + y/4 - z/2)*Hf_N2;
DH_comb = DH_comb - (Hf_CxHy + (x + y/4 - z/2)*(Hf_O2 + 3.78*Hf_N2));
PCI = -DH_comb/MM_carb;
psi_s = (x+y/4)*(2*MM_O+3.76*2*MM_N)/MM_carb;
Q_p = phi * PCI / (phi + psi_s);

%% Exportation des resultats

save('input_variables/proprietes_gaz.mat')

%% plots

% hold on
% plot(T,cv_air,'LineWidth',2)
% plot(T,cv_frais,'LineWidth',2)
% plot(T,cv_prod,'LineWidth',2)
% legend('Air', 'Frais', 'Brules','Location','Northwest') %Northwest / Northeast
% xlabel('Temperature (K)','fontweight','bold')
% ylabel('c_v (J/kg-K)','fontweight','bold')
% title('c_v=f(T)')
% % ylabel('\gamma_p','fontweight','bold')
% % title('\gamma_p=f(T)')
% grid()


% ft = fittype('a+b*x');
% f = fit(T(1:3)',cp_frais(1:3)',ft);
% plot(f)
% hold on
% plot(T(1:3),cp_frais(1:3),'*','LineWidth',2)
% xlabel('Temperature (K)','fontweight','bold')
% ylabel('c_p (J/kg-K)','fontweight','bold')
% title('c_p frais = f(T)')
% legend('Fitted', 'Points','Location','Northeast')
% grid()
