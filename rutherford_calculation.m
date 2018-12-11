% 1. RUTHERFORD INSTRUCTION MANUAL QUESTIONS
% divers
elementar = 1.6021766208*10^-19;
Zalpha = 2;
ZAu = 79;
chbar = 3.16152649*10^-26; %=c\cdot \hbar
epsilon2 = 1.44*10^-9 * 1.6021766208*10^-19; % in Vm
epsilon0 = 8.85418781762 *10^-12; 
finestructure = epsilon2 / chbar;
r0 = 1.3*10^-13*10^-3; % in m; equation (7.3)
A1 = 2 + 2; % nukleonenzahl von alpha particle
A2 = 79 + 79; % nukleonenzahl von nucleus
mAm = 241.056822944; % u; mass of Am particle
mAu = 196.9665695; % in u; mass of Au particle
mNp = 237.048167253; % u; mass of Np particle
mHe = 4.00263250; % u; mass of He particle
malpha = 4.001487900; % u; mass of \alpha particle
uc2 = 931.49432; % in MeV
K1 = malpha / mAu;

% a): calculate T_{alpha i}:
% use the formula and table on p.18 and data from p.19:
Q0 = (mAm - mNp - mHe); % in u; energy of decay / without error
Igiven = 0.34 + 0.22 + 84.5 + 13.0 + 1.6; % assumed to be without / negliable error
    % =99.6600 percent: normalize:
Inormed = [0.34, 0.22, 84.5, 13.0, 1.6]*100/Igiven; % negliable error
    %=[0.3412    0.2208   84.7883   13.0444    1.6055]
EAi = [0, 33.20, 59.54, 102.96, 158.51]*10^-3; % in MeV; excitation energy; given values in instruction manual p18
Talphai = zeros(1,5);
err_Talphai =  zeros(1,5); % Talphai does not depend on any measured quantity
for i = 1:5
    Talphai(i) =(Q0*uc2 - EAi(i)) / (1+ malpha / mNp);
end
Talphai; %=[5.5174    5.4848    5.4589    5.4162    5.3616] MeV

% b): Tm with given transition probabilities:
% use table on p.18:
Tm = 0;
err_Tm =0;
summe =0;
for i=1:5
    summe = Talphai(i) * Inormed(i)/100;
    Tm = Tm + summe; 
end
Tm; % = 5.4520 MeV

% c): Tmc mean kinetic energy of alpha particles in COM-system:
Tmc = Tm / (1+K1)^2; %=5.2371 MeV; without error since err_Tm =0 and K1 given
err_Tmc = 0;

% d): minimal distance alpha-nucleus using eq. (7.1):
theta = pi; % for min. distance: b=0 and so \theta = \pi
K2 = ((K1* cos(theta) + sqrt(1- K1 * sin(theta)^2)) / (1 + K1))^2; %kinematic factor; without error
E0 = Tm * 10^6 * elementar; % in V; without error
Dmin = Zalpha * ZAu * epsilon2 * (1 + K1) / E0; %=4.2579e-14 m; without error
    
% e): differential cross section area for b=0 using eq. (4.13):
sigma = 1.296 * (Zalpha * ZAu / Tm)^2 * (sin(theta /2)^-4 - 2*(malpha / mAu)^2);
    %=1.0875e+03 (mb)/(sr), without error

% f): scattering angle \Theta and impact parameter b for x_{min}, x_{max}:
xmin = 5*10^-2; %m; without error
xmax = 15*10^-2; %m; without error
Thetamin = 0.52; %rad
Thetamax = 0.83; % rad
bmax = Zalpha * ZAu * epsilon2 *cot(Thetamin /2) / 2E0; %=6.8515e-26 m
bmin = Zalpha * ZAu * epsilon2 *cot(Thetamax /2) / 2E0; %=4.1368e-26 m
impact = [bmin, bmax];

% g): deviations from Rutherford formula for used energy range:
% here we calculate the minimal energy E_r for which we expect a deviation
% both without error
Dmin_deviation = (r0 * (A1^(1/3) + A2^(1/3))); %= 9.0916e-16 m
Er = Zalpha * ZAu * epsilon2 / (r0 * (A1^(1/3) + A2^(1/3))); %=4.0095e-11 V

% h): stopping power; gold - using Bethe-Bloch formula 
%   in eV/(10^15 Atome/cm^2) and keV/mu m:
KAu = -1.037; % correction
EBAu = 1059.81; % in eV
BCAu_eV =  (3.80 / (Tm) * ZAu * log((548.58 * Tm) /(EBAu) - KAu));
    %= 74.3576
BCAu_keV = BCAu_eV *10^-3*(5.91*10^22) *10^-15 * 10^-4; 
    % =439.45 keV/mu m; should be: 522keV/mu m
    
% i): 

% j): E_min and E_max for x_min, x_max of alpha particles with respect of
% finite dimension of detector and foil
R1 = 23*10^-3; % m; inner radius of foil
R2 = 27*10^-3; % m; outer radius of foil
delta = 73*10^-3; % m; distance source-foil
x = [xmin, xmax]+ 3*10^-3; % distances between foil and detector
AD = 50*10^-6; %m^2; area of the detector
rD = sqrt(AD / pi); % m; readius of the detector
a = sqrt(delta^2 + (R1)^2); % =0.0765
b =  sqrt(delta^2 + (R2)^2); % =0.0778
a1 = acos(R1 / a); %=1.2656 rad
a3 =  acos(R2 / b); %=1.2165 rad
a2 = [0,0]; %=[1.1935    1.4344] rad
a4 = [0,0];%=[1.0702    1.3835] rad
for i=1:2
    a2(i) = atan((x(i))/(R1 - rD/2));
    a4(i) = atan((x(i))/(R2 + rD/2));
end
thetadown = a1 + a2; %=[min: 2.4590, max: 2.6999] rad
thetaup = a3 + a4; %=[min: 2.2867, max: 2.600] rad
Edown = [0,0];
Eup = [0,0];
for i=1:2
    Edown(i) = Zalpha * ZAu * epsilon2 / (4*impact(i)) *cot(thetadown(i)/2); 
    Eup(i) = Zalpha * ZAu * epsilon2 / (4*impact(i)) *cot(thetaup(i)/2); 
end
Edown; % = [max: 0.0782; min: 0.0299] MeV
Eup; % = [max: 0.1003; min: 0.0369] MeV

% k): stopping power: air - using Bethe-Bloch
Kair = 0.710; % correction
EBair = 94.22; % in eV
ZN = 7;
BCair_eV =  (3.80 / (Tm) * ZN * log((548.58 * Tm) /(EBair) - Kair));
    %= 74.3576
BCair_keV = BCair_eV *10^-3*(5.91*10^22) *10^-15 * 10^-4;
    % =439.45 keV/mu m; should be: 522keV/mu m

% l): correction for \sigma(\theta) with respect to finite mass of gold
% nucleus: we have to use formula (4.10) with E = Tmc: in the centre of
% mass system:
E_corr = Tmc * 10^6 * elementar; %V
sigma_min_corr = ((Zalpha * ZAu * epsilon2)/(4* E_corr^2))^2 * 1/sin(Thetamax / 2)^4;
    %=0.006340318283385
sigma_max_corr = ((Zalpha * ZAu * epsilon2)/(4* E_corr^2))^2 * 1/sin(Thetamin / 2)^4;
    %=0.038359494392539

% m) differential cross section area for xmin; with respect to finite
% dimension of foil and detector
% note: xmin <-> thetamax; and Tm [MeV]-> E [V]
E = Tm * 10^6 * elementar; %V
sigma_min = ((Zalpha * ZAu * epsilon2)/(4* E^2))^2 * 1/sin(Thetamax / 2)^4;
    %=0.0054 
    
  % correction is given by dividing (4.10)/ (4.12) from the instr. manual:
  corrmin = 1/(1-2*(malpha / mAu)^2 * sin(Thetamin / 2)^4);
    %=1.000003605507719
  corrmax = 1/(1-2*(malpha / mAu)^2 * sin(Thetamax / 2)^4);
    %=1.000021814042338

    
% 2. RUTHERFORD ERROR CALCULATION 
% activity of Am-241 source
delta = 73*10^-3; % m
Rhole = 0.5*10^-3; % m
err_Rhole = 0.1*10^-3; % m
t = 90.09; % s
err_t = 0.3; % s
ct = 3562; % counts
err_ct = 60; % ct
Ahole = pi * Rhole^2;
err_A_hole = pi * 2*Rhole *err_Rhole;
omega_tot = 4*pi;
omega_hole = Ahole / delta^2;
err_omega_hole = err_A_hole / delta^2;
    activity = ct / t * omega_tot / omega_hole; %=3.3712*10^6 Bq of Am-241 source
    err_activity = sqrt((omega_tot * err_ct / (t * omega_hole))^2 + ((t * omega_tot * err_t)/(t^2 * omega_hole))^2 );  % =5.6787e+04 Bq
            %+ ((t * omega_tot * err_omega_hole)/(t * omega_hole^2))^2:
            %negliable since err_omega_hole is small comparative to the other errors
    

% 3. RUTHERFORD CALCULATION / PLOTS
fname = 'C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Rutherford\plots';

% 3.1. DISCRIMINATOR CURVE
[~, ~, raw] = xlsread('C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Rutherford\data\discriminator_curve.xlsx','Tabelle1','B2:D17');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,3]);
raw = raw(:,2);
data = reshape([raw{:}],size(raw));
turns = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.75];
counts = cellVectors(:,1);
err_counts = data(:,1);
sec = cellVectors(:,2);
err_sec = 0.3;
clearvars data raw cellVectors;

counts=  str2double(counts).';
err_counts = err_counts.'; 
sec = str2double(sec).';
err_turns = 0.05 * ones(1, 16);
counts_sec =0;
err_counts_sec =0;
for i=1:16
    counts_sec(i) = counts(i) / sec(i);
    err_counts_sec(i) = sqrt((1/sec(i) * err_counts(i))^2 + (err_sec * counts(i)/sec(i)^2)^2);
end
counts_sec;
err_counts_sec;
y = erf(turns);
xinterval = linspace(-1, 6);
%yinterval = [0:0.01:40];

% create the erfc-fit:
[xData, yData] = prepareCurveData( turns, counts_sec );
ft = fittype( 'a*erf(b*x + c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-20 1.2 -3.5 20]; % startpoints for [a b c d]
[fitresult, gof] = fit( xData, yData, ft, opts );
Discr_coeff = coeffvalues(fitresult);
        %ft = fittype('a*sin(b*x)');
        %fitt1 = fit(xdata,ydata,ft,'StartPoint',[-20 1.2]);
        %[d1,d2] = differentiate(fitt1,xdata);

% create the corresponding gaussian :
% gaussian = (Discr_coeff(1) * Discr_coeff(2) * 2 / sqrt(pi)) .* exp(-(Discr_coeff(2) * xData(:) + Discr_coeff(3)).^2)
fitt1 = fit(xData, yData, ft, 'Startpoint', [-20 1.2 -3.5 20]);
[d1, d2] = differentiate(fitt1, xinterval);

% plot discriminator curve with x axis: turns:
figure( 'Name', 'erf' );
h = plot( fitresult, xData, yData );
hold on
errorbar(turns(:), counts_sec(:), err_counts_sec(:)/2, err_counts_sec(:)/2, err_turns(:)/2, err_turns(:)/2, 'bo')
hold on
% here the fit line:
plot(xinterval(:), -d1(:), 'g-');
grid on
legend('data', 'erfc', 'data', 'gaussian', 'Location', 'best' );
title('discriminator curve')
xlabel('turns')
ylabel('counts per second')
saveas(gcf, fullfile(fname, 'discriminator_curve.eps'), 'epsc');

% 4. ANGULAR DISTRIBUTION
% 4.1. Import the data:
    opts = spreadsheetImportOptions("NumVariables", 7);
    opts.Sheet = "Tabelle1";
    opts.DataRange = "A2:G13";
    opts.VariableNames = ["x_tildecm", "thetarad", "z_counts", "z_times", "z_U_counts", "z_U_times", "err_z_counts"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7], "EmptyFieldRule", "auto");
    tbl = readtable("C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Rutherford\data\Winkelverteilung.xlsx", opts, "UseExcel", false);
    x_tildecm = tbl.x_tildecm;
    thetarad = tbl.thetarad;
    z_counts = tbl.z_counts;
    z_times = tbl.z_times;
    z_U_counts = tbl.z_U_counts; % 'U' refers to 'underground'
    z_U_times = tbl.z_U_times;
    err_z_counts = tbl.err_z_counts; % is square root of z_counts
    clear opts tbl    
    
    % convert the imported data from string to double:
    x_tilde = str2double(x_tildecm); % m; 3mm are included
    x_tilde_m = x_tilde * 10^-2;
    theta_rad = str2double(thetarad); % rad
    z_counts = str2double(z_counts); % #
    z_times = str2double(z_times); % s
    z_U_counts = str2double(z_U_counts); % #
    z_U_times = str2double(z_U_times); % s
    err_z_counts = str2double(err_z_counts); % #; counting error
    
    err_z_times = 0.3 * ones(12, 1);
    err_z_U_times = 0.3 * ones(12, 1);
    err_x_tilde_m = 0.001 * ones(12,1);
    
 % 4.2. plot: sin(Theta/2) gegen Na / (OmegaD*t)
    % we have to determine Na (page 10):
    OmegaF = 0.0998; % given in equation (9.7)
    RA = (R1 + R2)/2; % mean radius of the gold foil
    OmegaD = zeros(12, 1); % =[0.1997, 0.1732, 0.1485, 0.1267, 0.1074, 0.0900, 0.0747, 0.0612, 0.0493, 0.0389, 0.0126, 0.0222]e-05
    err_OmegaD = zeros(12, 1);
    X = zeros(12, 1); % x axis in the following 2 plots
    err_X = zeros(12, 1);
    Y = zeros(12, 1); % y axis in the following 2 plots; note: Na / t corresponds to counts/second
    err_Y = zeros(12, 1);
    err_theta_rad = 0.01;
    Ct_s_angular = zeros(12, 1);
    err_z_counts = zeros(12, 1);
    err_ct_s_angular = zeros(12, 1);
    for i = 1:12
        Ct_s_angular(i) = z_counts(i) / z_times(i);
        err_z_counts(i) = sqrt(z_counts(i)); % counting error
        err_ct_s_angular(i) = sqrt((err_z_counts(i) / z_times(i))^2 + (z_counts(i) * err_z_times(i) / (z_times(i))^2)^2);    
        OmegaD(i) = pi * rD^2 * x_tilde_m(i) /((x_tilde_m(i))^2 + RA^2)^(3/2); % with error: x_tilde_m / RA and rd are given
        err_OmegaD(i) = sqrt(((pi * rD^2 * (RA^2-2*(x_tilde_m(i))^2))/(RA^2 + (x_tilde_m(i))^2)^(5/2))^2 * (err_x_tilde_m(i))^2);
        X(i) = sin(theta_rad(i) / 2); % with error: theta_rad
        err_X(i) = sqrt((cos(theta_rad(i) / 2 * err_theta_rad))^2);
        Y(i) = Ct_s_angular(i) / OmegaD(i); % with error: counts_sec, OmegaD
        err_Y(i) = sqrt((err_ct_s_angular(i) / OmegaD(i))^2 + (Ct_s_angular(i) / OmegaD(i)^2 * err_OmegaD(i))^2);
    end
    
    % CARTESIAN PLOT:
    xspace1 = linspace(0.24, 0.42) ;
        %for the fitline:
        [XXData, YYData] = prepareCurveData( X, Y );
        FT = fittype( 'a / x^4', 'independent', 'x', 'dependent', 'y' );
        OPTS = fitoptions( 'Method', 'NonlinearLeastSquares' );
        OPTS.Display = 'Off';
        OPTS.StartPoint = 0.132420604790519;
        [FITRESULT, GOF] = fit( XXData, YYData, FT, OPTS );
    figure( 'Name', 'erf' );
    h = plot( FITRESULT, XXData, YYData ); 
    hold on
    errorbar(X(:), Y(:), err_Y(:),  'bo')
    grid on
    xlabel('sin(\theta / 2)')
    ylabel('N_a / (\Omega_D t)')
    title('angular distribution')
    H = legend('data', 'fit', 'Location', 'Best');
    saveas(gcf, fullfile(fname, 'angular_distribution.eps'), 'epsc');
    
  % 4.3 plot: Theta vs. Na / (OmegaD*t); in the loglog-plot we expect a linear dependence (cf. p10)
    % LOGLOG - PLOT:
    XX = log(X(:)); % = ln(sin(theta / 2))
    YY = log(Y(:)); % = ln(N_a / (Omega_D * t))
    err_YY = sqrt((err_Y(:) ./ Y(:)).^2);
              
    [XX, YY] = prepareCurveData( XX, YY );
    fFt = fittype( 'a*x + b', 'independent', 'x', 'dependent', 'y' );
    oOpts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    oOpts.Display = 'Off';
    oOpts.StartPoint = [0.0579159481171146 0.304881922129305];
    [fFitresult, gGof] = fit( XX, YY, fFt, oOpts );
    coeffvals = coeffvalues(fFitresult); %=[-4.3399   -0.2400] = [a, b]
    err_b = (0.3345 + 0.1455) /2; %positive, negative direction
    err_a = (0.2710 + 0.2969)/2; %positive, negative direction
       
    xspace2 = linspace(-1.4, -0.9) ;
    figure
    errorbar(XX(:), YY(:), err_YY(:), 'bo')
    hold on
    grid on
    BBB = polyval(coeffvals, xspace2);
    plot(xspace2, BBB, 'c-')
    ylim([3.5 6])
    xtext = -1.208;
    ytext = 5.003; 
    str = {'\leftarrow y =-4.34x -0.24'};
    text(xtext, ytext, str);
    xlabel('ln(sin(\theta / 2))')
    ylabel('ln(N_a / (\Omega_D t))')
    G = legend('data', 'linear fit', 'Location', 'Best');
    title('angular distribution (loglog scale)')
    saveas(gcf, fullfile(fname, 'angular_distribution_loglog.eps'), 'epsc');
        
    C = exp(coeffvals(2)); %= 0.7866
    epsilon2_tilde = 1.44;
    foil_thickness = 1*10^-6; %m, Skript: 5.2.2.
    n_AK = (5.91*10^22 * 10^6) * (mAu * 10^-3) / (6.02214085774 * 10^23); %=19330 kg/m^3; literature: 19320kg/m^3 [4]
    estimate_activity = (4 *pi * C)/(n_AK* foil_thickness * OmegaF) * ((Zalpha * ZAu * epsilon2)/(4 * Tm  * elementar *10^6))^-2; 
        %= 4.7078e+07 Bq according to eq (4.13): Tm in MeV
    % error calculation:
    err_C = sqrt((exp(b) * err_b)^2); % =0.2594
    err_estimate_activity = sqrt((4*pi*err_C /(OmegaF *foil_thickness* n_AK) * (Zalpha*ZAu*epsilon2_tilde / (4*Tm*10^6))^-2)^2 + ...
        (((4*pi*C)/(OmegaF*foil_thickness * n_AK))*2*((4*Tm*10^6)/(Zalpha*ZAu * epsilon2_tilde))*((4*err_Tm)/(Zalpha*ZAu*epsilon2_tilde)))^2);
        % = 1.5526e+07 Bq
   
    
% 5. CHI SQAURE TEST
h_angular = chi2gof(Y); %=0
    % h = chi2gof(x) returns a test decision for the null hypothesis that the data in vector x comes from a normal 
    % distribution with a mean and variance estimated from x, using the chi-square goodness-of-fit test. 
    % The alternative hypothesis is that the data does not come from such a distribution. 
    % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
        % source: https://ch.mathworks.com/help/stats/chi2gof.html
