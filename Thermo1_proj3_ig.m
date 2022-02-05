clear all, close all, clc

global R A B C D;
R = 8.314462;   % l*kPa / mol*K or J/mol*K
A = 5.457; B = 1.045 * 10^-3; C = 0; D = -1.157 * 10^5;

P_initial = 101.325 * 1;  % kPa
P_final = 11000;          % kPa

T1 = 25 + 273.15; % ambient temperature
Tref = -30 + 273.15; % refrigerated for storage (low pressure)

EperCO2 = 10.5*0.315; % MJ/kg

fprintf(' << Calculation assuming ideal-gas state >>\n\n')
fprintf(' Initial temperature        %.2f K\n', T1)
fprintf(' Initial pressure           %.3f kPa\n', P_initial);
fprintf(' Final pressure             %.1f kPa\n', P_final);

n = 3;  % number of compressors in the chain
r_total = P_final / P_initial ; % total compression ratio
r = r_total^(1/n);              % single compression ratio
fprintf(' Number of compressors      %d\n', n);
fprintf(' Total compression ratio    %.3f\n', r_total);
fprintf(' Single compression ratio   %.3f\n', r);

eff_comp = 0.65;    % single compressor efficiency
fprintf(' Compressor efficiency      %.f %% \n', 100*eff_comp);

P = zeros(n+1,1);
for i = 1:n+1
    P(i) = P_initial * r^(i-1);
end

fprintf(' Ambient temperature        %.2f K\n', T1);
fprintf(' Refrigeration temperature  %.2f K\n', Tref);

W_comp_total = 0;
W_other = 0.970*44.01*1000;

fprintf('\n\n >> Compression to %d MPa\n', P_final/1000);

for j = 1:n
    
    T1 = 25 + 273.15;

    
    P1 = P(j);
    P2 = P(j+1);

    T2_isen = zeros(100,1);
    T2_isen(1) = T1 + 100;
    for i = 1:10000
        fprintf('\nT2(isen)[%d] = %.3f',i, T2_isen(i));
        T2_isen(i+1) = T1*(P2/P1)^(1/MCPS(T1,T2_isen(i)));
        if abs( (T2_isen(i+1) - T2_isen(i)) / T2_isen(i) ) < 0.0001
            T2_is = T2_isen(i+1);
            clear T2_isen
            break
        end
    end
        fprintf('\nT2(isen)[%d] = %.3f',i+1, T2_is);
        
    fprintf('\n');
    CpH = R * MCPH(T1,T2_is);   % heat capacity
    Ws_comp_isen = CpH * (T2_is - T1);  % minimum work of a compressor
    Ws_comp = Ws_comp_isen / eff_comp;  % actual required work of a compressor

    T2_real = zeros(100,1);
    T2_real(1) = T2_is;
    for i = 1:10000
        
        T2_real(i+1) = T1 + Ws_comp / (R*MCPH(T1,T2_real(i)));
        fprintf('\nT2(real)[%d] = %.3f',i, T2_real(i));
        if abs( (T2_real(i+1) - T2_real(i))/T2_real) < 0.0001
            T2_rl = T2_real(i+1);
            clear T2_real
            break
        end
    end
        fprintf('\nT2(real)[%d] = %.3f',i+1, T2_rl);

    W_comp_total = W_comp_total + Ws_comp;
    fprintf('\n\n - Compression step %d: ',j)
    fprintf('%.3f kPa to %.3f kPa', P1, P2);
end

fprintf('\n');
fprintf('\n       Isentropic work  %.1f J/mol', n*Ws_comp_isen);
fprintf('\n Real compression work  %.1f J/mol', W_comp_total);
fprintf('\n                        %.5f MJ/kg', W_comp_total/44.01/1000);
fprintf('\n             Highest T  %.3f K\n', T2_rl);

if P_final < 7000
    
    W_ref = 0.1264*44.01*1000; % molar mass of CO2: 44.01 g/mol
    fprintf('\n >> Refrigeration\n');
    fprintf('\n                 T_hot  %.3f K',T1);
    fprintf('\n                T_cold  %.3f K', Tref);
    fprintf('\n    Refrigeration work  %.2f J/mol', W_ref);
    fprintf('\n                        %.5f MJ/kg', W_ref/44.01/1000);
    fprintf('\n        Coefficient of')
    fprintf('\n           Performance  %.2f\n', Tref/(T1-Tref))  
else
    fprintf('\n');
    W_ref = 0;
end

fprintf('\n >> Total work required\n');
fprintf('\n        Compression %.1f J/mol',W_comp_total);

if P_final<7000
fprintf('\n      Refrigeration  %.1f J/mol',W_ref);
end

fprintf('\n              Other %.1f J/mol',W_other);
fprintf('\n              Total %.1f J/mol',W_other+W_comp_total+W_ref);
fprintf('\n                    %.4f MJ/kg',(W_other+W_comp_total+W_ref)/44.01/1000);
fprintf('\n Emission intensity %.4f MJ/kg',EperCO2);

fprintf('\n\n >> Energy penalty\n');
fprintf('\n %.2f %%',(W_other+W_comp_total+W_ref)/44.01/1000/EperCO2*100);
fprintf('\n\n');

function MCPS = MCPS(T1,T2)
    global A B C D;
    MCPS = A + (B + ( C + D/(T1^2*T2^2) )*((T1+T2)/2))*(T2-T1)/log(T2/T1);
end

function MCPH = MCPH(T1,T2)
    global A B C D;
    MCPH = A + (B/2)*(T2+T1) + (C/3)*(T2^2 + T1^2 + T1*T2) + D/(T1*T2);
end

