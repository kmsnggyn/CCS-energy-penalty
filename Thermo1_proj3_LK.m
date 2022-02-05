clear all, close all, clc

global R AA BB CC DD;
R = 8.314462;   % l*kPa / mol*K or J/mol*K
AA = 5.457; BB = 1.045 * 10^-3; CC = 0; DD = -1.157 * 10^5; %Cp_ig const. for CO2

P_initial = 101.325;   % kPa
P_final = 11000;       % kPa

T1 = 25 + 273.15; % ambient temperature
Tref = -30 + 273.15; % refrigerated for storage (low pressure)

EperCO2 = 10.5*0.315; % MJ/kg

fprintf(' << Calculation using Lee-Kelser EOS >>\n\n')
fprintf(' Initial temperature        %.2f K\n', T1)
fprintf(' Initial pressure           %.3f kPa\n', P_initial)
fprintf(' Final pressure             %.2f kPa\n', P_final);

n = 3;
r_total = P_final / P_initial ;
r = r_total^(1/n);
fprintf(' Number of steps            %d \n', n)
fprintf(' Total compression ratio    %.3f \n', r_total)
fprintf(' Single compression ratio   %.3f \n', r)

eff_comp = 0.65;
fprintf(' Compressor efficiency      %.f %% \n', 100*eff_comp);

P_comp = zeros(n+1,1);
for i = 1:n+1
    P_comp(i) = P_initial * r^(i-1);
end

fprintf(' Ambient temperature        %.2f K\n', T1);
fprintf(' Refrigeration temperature  %.2f K\n', Tref);

W_other = 0.970*44.01*1000;

fprintf('\n\n >> Compression to %d MPa\n', P_final/1000);

w = 0.22394; w0 = 0.3978; Tc = 304.1282; Pc = 7377.3;
T2(1) = 500;
HR = zeros(2,1);
SR = zeros(2,1);
a = 1;

for j = 1:1000
    Tmat = [T1, T2(j)];
    Pmat = [P_comp(1), P_comp(2)];

    for i=1:2 % To calculate with two sets of T and P
        T = Tmat(i);
        P = Pmat(i);
        Tr = T/Tc;
        Pr = P/Pc;

        % Calculating for simple fluids
        b1 = 0.1181193;
        b2 = 0.265728;
        b3 = 0.154790;
        b4 = 0.030323;
        c1 = 0.0236744;
        c2 = 0.0186984;
        c3 = 0.0;
        c4 = 0.042724;
        d1 = 0.155488*10^-4;
        d2 = 0.623689*10^-4;
        beta = 0.65392;
        gamma = 0.060167;
        B = b1 - b2/Tr - b3/Tr^2 - b4/Tr^3;
        C = c1 - c2/Tr + c3/Tr^3;
        D = d1 + d2/Tr;

        % Numerically solving the EOS to obtain Vr (molar volume)
        syms Vr
        LK_eos = Pr*Vr/Tr == 1 + (B/Vr) + (C/Vr^2) + (D/Vr^5) + (c4/(Tr^3*Vr^2))*(beta + gamma/Vr^2)*exp(- gamma/Vr^2);
        Vr0 = vpasolve(LK_eos,Vr);
        Vr = Vr0;

        % Calculating the compressibility for simple fluids
        Z0 = Pr*Vr0/Tr;
        Z = Z0;

        % Calculating the residual properties for simple fluids
        E = (c4/(2*Tr^3*gamma))*(beta+1-(beta+1+gamma/Vr^2)*exp(-gamma/Vr^2));
        HR_RTc_0 = Tr*(Z-1-(b2 + 2*b3/Tr + 3*b4/Tr^2)/(Tr*Vr) - (c2- 3*c3/Tr^2)/(2*Tr*Vr^2) + d2/(5*Tr*Vr^5)+3*E);
        SR_R_0 = log(Z) - ((b1 + b3/Tr^2 + 2*b4/Tr^3 )/Vr) - ((c1- (2*c3)/(Tr^3))/(2*Vr^2)) - (d1/(5*Vr^5)) + 2*E;
        
        % Calculations for reference fluids
        b1 = 0.2026579;
        b2 = 0.331511;
        b3 = 0.027655;
        b4 = 0.203488;
        c1 = 0.0313385;
        c2 = 0.0503618;
        c3 = 0.016901;
        c4 = 0.041577;
        d1 = 0.48736*10^-4;
        d2 = 0.0740336*10^-4;
        beta = 1.226;
        gamma = 0.03754;

        B = b1 - b2/Tr - b3/Tr^2 - b4/Tr^3;
        C = c1 - c2/Tr + c3/Tr^3;
        D = d1 + d2/Tr;

        % Numerically solving the EOS to obtain Vr (molar volume)
        syms Vr
        LK_eos = Pr*Vr/Tr == 1 + (B/Vr) + (C/Vr^2) + (D/Vr^5) + (c4/(Tr^3*Vr^2))*(beta + gamma/Vr^2)*exp(-gamma/Vr^2);
        Vr1 = vpasolve(LK_eos,Vr,10);
        Vr = Vr1;

        % Calculating the compressibility factor for reference fluids
        Z1 = Pr*Vr1/Tr;
        Z = Z1;

        % Calculating the residual properties for reference fluids
        E = (c4/(2*Tr^3*gamma))*(beta+1-(beta+1+gamma/Vr^2)*exp(-gamma/Vr^2));
        HR_RTc_1 = Tr*(Z-1-(b2 + 2*b3/Tr + 3*b4/Tr^2)/(Tr*Vr) - (c2- 3*c3/Tr^2)/(2*Tr*Vr^2) + d2/(5*Tr*Vr^5)+3*E);
        SR_R_1 = log(Z) - ((b1 + b3/Tr^2 + 2*b4/Tr^3 )/Vr) - ((c1-(2*c3)/(Tr^3))/(2*Vr^2)) - (d1/(5*Vr^5)) + 2*E;

        % Calculating the properties for the substance of interest
        Z = Z0 + (w/w0)*(Z1-Z0);
        Vr = Z*Tr/Pr;

        HR(i) = R*Tc* ( HR_RTc_0 + (w/w0) * (HR_RTc_1 - HR_RTc_0) );
        SR(i) = R* ( SR_R_0 + (w/w0) * (SR_R_1 - SR_R_0) );

    end

    delta_H = - HR(1) + HR(2) + R * MCPH(Tmat(1),Tmat(2)) * (Tmat(2)-Tmat(1));
    delta_S = - SR(1) + SR(2) + R * MCPS(Tmat(1),Tmat(2)) * log(Tmat(2)/Tmat(1)) - R*log(r);

    S_tot(j) = delta_S;
%     fprintf('.')

    fprintf('\n T2[%d] = %.4f, delta S = %.4f ', j, T2(j), S_tot(j));
    
    if j > 1 && S_tot(j) * S_tot(j-1) < 0
        a = a+1;
        if a > 6
            break
        end
    end

    if S_tot(j) > 0
        T2(j+1) = T2(j) - 10^(2-a);
    elseif delta_S < 0
        T2(j+1) = T2(j) + 10^(2-a);
    end
end

W_comp_isen = n * delta_H;
W_comp = W_comp_isen / eff_comp;

S_tot2 = n * delta_S;

fprintf('\n');
fprintf('\n           Change in S  %.4f J/mol K', S_tot2);
fprintf('\n           Change in H  %.1f J/mol \n', W_comp_isen);
fprintf('\n       Isentropic work  %.1f J/mol', W_comp_isen);
fprintf('\n Real compression work  %.1f J/mol', W_comp);
fprintf('\n                        %.5f MJ/kg', W_comp/44.01/1000);
fprintf('\n             Highest T  %.3f K\n', T2(j));

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
    W_ref = 0;
    fprintf('\n');
end

fprintf('\n >> Total work required\n');
fprintf('\n        Compression %.1f J/mol',W_comp);

if P_final<7000
fprintf('\n      Refrigeration  %.1f J/mol',W_ref);
end

fprintf('\n              Other %.1f J/mol',W_other);
fprintf('\n              Total %.1f J/mol',W_other+W_comp+W_ref);
fprintf('\n                    %.4f MJ/kg',(W_other+W_comp+W_ref)/44.01/1000);
fprintf('\n Emission intensity %.4f MJ/kg',EperCO2);

fprintf('\n\n >> Energy penalty\n');
fprintf('\n %.2f %%',(W_other+W_comp+W_ref)/44.01/1000/EperCO2*100);
fprintf('\n\n');


function MCPS = MCPS(T1,T2)
    global AA BB CC DD;
    MCPS = AA + (BB + ( CC + DD/(T1^2*T2^2) )*((T1+T2)/2))*(T2-T1)/log(T2/T1);
end

function MCPH = MCPH(T1,T2)
    global AA BB CC DD;
    MCPH = AA + (BB/2)*(T2+T1) + (CC/3)*(T2^2 + T1^2 + T1*T2) + DD/(T1*T2);
end

