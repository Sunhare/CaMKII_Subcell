%% Parameters for external modules

% ECC and CaM modules
freq = 1.0;                 % [Hz] CHANGE DEPENDING ON FREQUENCY
cycleLength = 1e3/freq;     % [ms]
CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM] 
CaNtotDyad = 3e-3/8.293e-4; % [uM] 
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

% ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
expression = 'WT';
CKIIOE = 0; % Should be zero during 'WT' and 'KO' runs

if strcmp(expression,'OE')
    CKIIOE = 1; % Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    CaMKIItotDyad = 120*6;          % [uM] 
    CaMKIItotSL = 120*8.293e-4*6;   % [uM]
    CaMKIItotCyt = 120*8.293e-4*6;  % [uM]
elseif strcmp(expression,'KO')
    CaMKIItotDyad = 0;              % [uM] 
    CaMKIItotSL = 0;                % [uM]
    CaMKIItotCyt = 0;               % [uM]
end

% Parameters for CaMKII module
LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
Ligtot = 0%; 0.1 or 0.02    % [uM] - SET LIGAND CONCENTRATION HERE
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = PLBtot;          % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]
PLMtotBA = 48;              % [uM] - [umol/L cytosol]
MyototBA = 70;              % [uM] - [umol/L cytosol]
IKrtotBA = 0.025;           % [uM] - [umol/L cytosol]
IClCatotBA = 0.025;         % [uM] - [umol/L cytosol]

% For Recovery from inactivation of LCC
recoveryTime = 10;  % initialize to smallest value
%% Collect all parameters and define mass matrix for BAR module

p = [cycleLength,recoveryTime,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,...
    PLMtotBA,MyototBA,IKrtotBA,IClCatotBA,CKIIOE];

% Need to define Mass matrix for BAR portion
m0 = ones(1,83+6+45+6);   % All state variables in other modules, not BAR
m1 =[0,     0,      0,      1,      1,      1,       1,       1,       1];
m2 =[0,     0,      0,      1,      1,      1];
m3 =[0,     0,      0];
m4 =[1,     1,      0,      0,      1,      1];
m5 =[1,     1,      0,      0,      1,      1,       1,       1      0,      0,      1,      1];

M = diag([m0,m1,m2,m3,m4,m5]); 
%% Establish and define globals

global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store 
global gates_store Jserca_store IKs_store Jleak_store ICFTR_store Incx_store
global I_NaK_store I_kr_store I_Nabk_store
global Lmyo_store Fmyo_store
global IKp_store I_clca_store I_clbk_store Ipca_store Icabk_store Cai_store

tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6); I_to_store=zeros(3,1e6);
I_Na_store = zeros(1,1e6); I_K1_store = zeros(1,1e6); ibar_store=zeros(1,1e6);
gates_store = zeros(2,1e6); Jserca_store = zeros(1,1e6); IKs_store = zeros(1,1e6);
Jleak_store = zeros(1e6,2); ICFTR_store = zeros(1,1e6); Incx_store = zeros(1,1e6);
I_NaK_store = zeros(1,1e6); I_kr_store = zeros(1,1e6); I_Nabk_store = zeros(1,1e6);

Fmyo_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6);

IKp_store = zeros(1,1e6); I_clca_store=zeros(1,1e6); I_clbk_store=zeros(1,1e6);
Ipca_store=zeros(1,1e6); Icabk_store=zeros(1,1e6); Cai_store=zeros(1,1e6);
%% Assign initial conditions

% Isometric contraction (mechFlag = 0 in ECC module)
load yfinal_myo_WT_1Hz_isometric                   % control (steady-state, 1-Hz)
%load yfinal_myo_WT_1Hz_isometric_240iso            % 100 nM [ISO] (240-s, 1-Hz) 
%load yfinal_myo_WT_1Hz_isometric_240iso_XBCa       % ISO-XBCa
%load yfinal_myo_WT_1Hz_isometric_240iso_XBcy       % ISO-XBcy
%load yfinal_myo_WT_1Hz_isometric_240iso_XBCa_XBcy  % ISO-XBCa-XBcy
%load yfinal_myo_WT_1Hz_isometric_240iso_titin      % ISO-titin
%load yfinal_myo_WT_1Hz_isometric_240iso_cytofl     % ISO-cytofl

% Isotonic contraction (mechFlag = 1 in ECC module)
% load yfinal_myo_WT_1Hz_isotonic                    % control (steady-state, 1-Hz)
%load yfinal_myo_WT_1Hz_isotonic_240iso             % 100 nM [ISO] (240-s, 1-Hz) 
%load yfinal_myo_WT_1Hz_isotonic_240iso_XBCa        % ISO-XBCa
%load yfinal_myo_WT_1Hz_isotonic_240iso_XBcy        % ISO-XBcy
%load yfinal_myo_WT_1Hz_isotonic_240iso_XBCa_XBcy   % ISO-XBCa-XBcy
%load yfinal_myo_WT_1Hz_isotonic_240iso_titin       % ISO-titin
%load yfinal_myo_WT_1Hz_isotonic_240iso_cytofl      % ISO-cytofl

y0n=yfinal;
%% Run single simulation

tspan = [0 3*1e3]; % [ms]
options = odeset('Mass',M,'RelTol',1e-5,'MaxStep',2); 

% % % % % % % % MAIN ODE % % % % % % % %
% [t,y] = ode15s(@rabbit_myofilament_masterODEfile,tspan,y0n,options,p);
% % % % % % % % MAIN ODE % % % % % % % %

% function dydt = rabbit_myofilament_masterODEfile(t,y,p)
% Negroni et al model
% This function calls the ode files for EC coupling, CaM reactions, CaMKII
% phosphorylation module, and PKA phosphorylation module

%% Collect params and ICs for each module

% Allocate ICs for each moduel
% y(60) -> y(65) are state transitions for mode 1 junctional LCCs
% y(66) -> y(71) are state transitoins for mode 2 junctional LCCs
% y(72) -> y(77) are state transitions for mode 1 sarcolemmal LCCs
% y(78) -> y(83) are state transitions for mode 2 sarcolemmal LCCs

% y(84) -> y(89) are state transitions for myofilament model

y_ecc = y(1:83+6);  % Ca_j is y(36), Ca_sl is y(37), Ca_cytosol is y(38) , Shannon Model, Restreppo in our case
y_camDyad = y(83+6+1:83+6+15); %14 ODEs 
y_camSL = y(83+6+15+1:83+6+15+15); %14 ODEs
y_camCyt = y(83+6+15+15+1:83+6+15+15+15); %14 ODEs
y_CaMKII = y(83+6+45+1:83+6+45+6);      % 6 state vars in CaMKII module
y_BAR = y(83+6+45+6+1:83+6+45+6+30+6);    % 30+6 state vars in BAR module

% break up parameters from master function into modules
paramsCell=mat2cell(p,ones(size(p,1),1),ones(size(p,2),1));
[cycleLength,recoveryTime,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot...
    PLMtotBA,MyototBA,IKrtotBA,IClCatotBA,CKIIOE]=paramsCell{:};

K = 135; % [mM]
Mg = 1;  % [mM]
%% Distribute parameters by module

% CaM module
CaDyad = y(36)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compart_dyad = 2;
% ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad];
CaSL = y(37)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartSL = 1;
pCaMSL = [K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength, compartSL];
CaCyt = y(38)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartCyt = 0;
pCaMCyt = [K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength, compartCyt];

% CaMKII phosphorylation module 
CaMKIIact_Dyad = CaMKIItotDyad.*(y(83+6+8)+y(83+6+9)+y(83+6+10)+y(83+6+11)); % Multiply total by fraction
CaMKIIact_SL = CaMKIItotSL.*(y(83+6+15+8)+y(83+6+15+9)+y(83+6+15+10)+y(83+6+15+11));
PP1_PLB_avail = y(83+6+45+6+22)./PP1_PLBtot + .0091;  % Active PP1 near PLB / total PP1 conc + basal value
pCaMKII = [CaMKIIact_Dyad,LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,...
           CaMKIIact_SL,LCCtotSL,PP1_SL,...
           PP1_PLB_avail];

LCC_CKdyadp = y(83+6+45+2)./LCCtotDyad;   %136 fractional CaMKII-dependent LCC dyad phosphorylation
RyR_CKp = y(83+6+45+4)./RyRtot;           %138 fractional CaMKII-dependent RyR phosphorylation
PLB_CKp = y(83+6+45+5)./PLBtot;           %139 fractional CaMKII-dependent PLB phosphorylation

% BAR (PKA phosphorylation) module
pBAR = [Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,PLMtotBA,MyototBA,IKrtotBA,IClCatotBA];
LCCa_PKAp = y(83+6+45+6+23)./LCCtotBA;
LCCb_PKAp = y(83+6+45+6+24)./LCCtotBA;
PLB_PKAn = (PLBtotBA - y(83+6+45+6+19))./PLBtotBA; % non-phosphorylated PLB targets
RyR_PKAp = y(83+6+45+6+25)./RyRtotBA;
TnI_PKAp = y(83+6+45+6+26)./TnItotBA;
IKs_PKAp = y(83+6+45+6+29)./IKstotBA;
ICFTR_PKAp = y(83+6+45+6+30)./ICFTRtotBA;
PLM_PKAp = y(83+6+45+6+31)./PLMtotBA;
Myo_PKAp = y(83+6+45+6+32)./MyototBA;
IKr_PKAp = y(83+6+45+6+35)./IKrtotBA;
IClCa_PKAp = y(83+6+45+6+35)./IClCatotBA;

% ECC module
pECC = [cycleLength,LCC_CKdyadp,RyR_CKp,PLB_CKp,...
        LCCa_PKAp,LCCb_PKAp,PLB_PKAn,RyR_PKAp,TnI_PKAp,IKs_PKAp,ICFTR_PKAp,...
        PLM_PKAp,Myo_PKAp,IKr_PKAp,IClCa_PKAp,CKIIOE,recoveryTime];
%% Solve dydt in each module

dydt_ecc = rabbit_myofilament_eccODEfile(t,y_ecc,pECC);
dydt_camDyad = rabbit_myofilament_camODEfile(t,y_camDyad,pCaMDyad);
dydt_camSL = rabbit_myofilament_camODEfile(t,y_camSL,pCaMSL);
dydt_camCyt = rabbit_myofilament_camODEfile(t,y_camCyt,pCaMCyt);
dydt_CaMKIIDyad = rabbit_myofilament_camkiiODEfile(t,y_CaMKII,pCaMKII);
dydt_BAR = rabbit_myofilament_barODEfile(t,y_BAR,pBAR);

% incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
global JCaCyt JCaSL JCaDyad;
dydt_ecc(36) = dydt_ecc(36) + 1e-3*JCaDyad;
dydt_ecc(37) = dydt_ecc(37) + 1e-3*JCaSL;
dydt_ecc(38) = dydt_ecc(38) + 1e-3*JCaCyt; 

% incorporate CaM diffusion between compartments
Vmyo = 2.1454e-11;      % [L]
Vdyad = 1.7790e-014;    % [L]
VSL = 6.6013e-013;      % [L]
kDyadSL = 3.6363e-16;	% [L/msec]
kSLmyo = 8.587e-15;     % [L/msec]
k0Boff = 0.0014;        % [s^-1] 
k0Bon = k0Boff/0.2;     % [uM^-1 s^-1] kon = koff/Kd
k2Boff = k0Boff/100;    % [s^-1] 
k2Bon = k0Bon;          % [uM^-1 s^-1]
k4Boff = k2Boff;        % [s^-1]
k4Bon = k0Bon;          % [uM^-1 s^-1]
CaMtotDyad = sum(y_camDyad(1:6))+CaMKIItotDyad*sum(y_camDyad(7:10))+sum(y_camDyad(13:15));
Bdyad = BtotDyad - CaMtotDyad; % [uM dyad]
J_cam_dyadSL = 1e-3*(k0Boff*y_camDyad(1) - k0Bon*Bdyad*y_camSL(1)); % [uM/msec dyad]
J_ca2cam_dyadSL = 1e-3*(k2Boff*y_camDyad(2) - k2Bon*Bdyad*y_camSL(2)); % [uM/msec dyad]
J_ca4cam_dyadSL = 1e-3*(k2Boff*y_camDyad(3) - k4Bon*Bdyad*y_camSL(3)); % [uM/msec dyad]
J_cam_SLmyo = kSLmyo*(y_camSL(1)-y_camCyt(1)); % [umol/msec]
J_ca2cam_SLmyo = kSLmyo*(y_camSL(2)-y_camCyt(2)); % [umol/msec]
J_ca4cam_SLmyo = kSLmyo*(y_camSL(3)-y_camCyt(3)); % [umol/msec]
dydt_camDyad(1) = dydt_camDyad(1) - J_cam_dyadSL;
dydt_camDyad(2) = dydt_camDyad(2) - J_ca2cam_dyadSL;
dydt_camDyad(3) = dydt_camDyad(3) - J_ca4cam_dyadSL;
dydt_camSL(1) = dydt_camSL(1) + J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL;
dydt_camSL(2) = dydt_camSL(2) + J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL;
dydt_camSL(3) = dydt_camSL(3) + J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL;
dydt_camCyt(1) = dydt_camCyt(1) + J_cam_SLmyo/Vmyo;
dydt_camCyt(2) = dydt_camCyt(2) + J_ca2cam_SLmyo/Vmyo;
dydt_camCyt(3) = dydt_camCyt(3) + J_ca4cam_SLmyo/Vmyo;
%% Collect all dydt terms

dydt = [dydt_ecc; dydt_camDyad; dydt_camSL; dydt_camCyt; dydt_CaMKIIDyad; dydt_BAR];
yfinal = y(end,:)';








