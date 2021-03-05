function [Final_cost_sep,Final_cost_collector,Final_cost_seals,Final_cost_pump,Final_cost_pcs,Final_cost_transformer,Final_cost_interconnection,Final_number_of_cells,Final_channel_height,Final_ASR,exitflag] = ssfbcostanalysis_0830_1

load ('Catholyte 4p2CB 11MnO2.mat');
load ('Cost components for power_future_upper');
load ('Reciprocating pump parameters');

n_cell_pump = 100;

%% Channel dimensions
channel_width = 2e-2; % 2cm

%% System requirement
power_output = 1000e3; %in W

%% Separator
ASR_sep = 0.1e-4; %R_sep 0.1 ohm.cm2 Hopkins et al.


%% Calculation

% Area specific impedance (in ohm.m2)
ASR_CT_pos = @(channel_height,n) r_ct_pos / (vol_frac_am_pos * channel_height * a_pos);
% ASR_CT_neg = @(channel_height,n) r_ct_neg / (vol_frac_am_neg * channel_height * a_neg);
ASR_CT_neg = 0;
ASR_CONT_pos = @(channel_height,n) r_cont_pos / (vol_frac_am_pos * a_pos);
ASR_CONT_neg = 0;
% ASR_CONT_neg = @(channel_height,n) r_cont_neg / (vol_frac_am_neg * a_neg);

ASR_OHM_pos = @(channel_height,n) (channel_height/2)*(1/catholyte_econd + 1/(catholyte_icond*vol_frac_elec_neg^1.5));
ASR_OHM_neg = @(channel_height,n) (channel_height/2)*(1/anolyte_econd + 1/(anolyte_icond*vol_frac_elec_neg^1.5));

ASR =@(channel_height,n) ASR_sep...
    + ASR_CT_pos(channel_height,n) + ASR_CT_neg...
    + ASR_CONT_pos(channel_height,n) + ASR_CONT_neg...
    + ASR_OHM_pos(channel_height,n) + ASR_OHM_neg(channel_height,n);


% Calculating A_s and delta_pressure_catholyte and anolyte
P_d =@(channel_height,n) (catholyte_phi_eq^2) * (1-eta_v) * eta_v / ASR(channel_height,n);   %Power density required for 1 cell; ASI is an input
A_s =@(channel_height,n) power_output / P_d(channel_height,n); %total area including multiple cells m2
L =@(channel_height,n) (A_s(channel_height,n) / (channel_width * n));
delta_pressure_catholyte =@(channel_height,n) 2 * (catholyte_tau_y/channel_height) * L(channel_height,n); % pressure drop for one cell
% delta_pressure_anolyte =@(vol_frac_am_pos,channel_height,n) 2 * (anolyte_tau_y/channel_height) * L(vol_frac_am_pos,channel_height,n); % pressure drop for one cell

% Calculating flow rate
flow_rate =@(channel_height,n) n_cell_pump * power_output / (n * catholyte_rho_e);

% Calculating efficiency
eps_pump =@(channel_height,n) 2 * (n/n_cell_pump) * delta_pressure_catholyte(channel_height,n) * flow_rate(channel_height,n) / (power_output * pump_eff);
eta_sys =@(channel_height,n) (1-eps_pump(channel_height,n));

% Calculating other costs
C_stack =@(channel_height,n) A_s(channel_height,n)*(Cost_sep + Cost_collector + Cost_seals) ; % in $/W
C_stack_n =@(channel_height,n) (C_stack(channel_height,n) / (power_output*eta_sys(channel_height,n))) * 1000; % in $/kW

pump_power_inHP = @(channel_height,n) delta_pressure_catholyte(channel_height,n) * flow_rate(channel_height,n) * 0.00134102; % in horsepower HP
Cost_pump =@(channel_height,n) (2753 * pump_power_inHP(channel_height,n)^0.57) * 1.59 * Pump_price_factor; % cost per pump
C_bop_n_pump =@(channel_height,n) (2*(n/n_cell_pump)*Cost_pump(channel_height,n) / power_output) * 1000; % in $/kW


% Calculating cost that scales with power
C_power =@(channel_height,n) Cost_BOP_no_pump + Cost_add + C_bop_n_pump(channel_height,n) + C_stack_n(channel_height,n) ; % in kW


lb = [2.5e-4 1000];
ub = [4e-3 14000];

IntCon = 2;


[sol,fval,exitflag,output] = ga(@(x) C_power(x(1),x(2)), 2, [], [], [], [], lb, ub, @(x) nonlinconstr_SSFB_0830_1(x,pump_power_inHP,pump_power_min_inHP), IntCon);

Final_cost_sep(1,1) = (Cost_sep/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2)); % in kW
Final_cost_collector(1,1) = (Cost_collector/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2))
Final_cost_seals(1,1) = (Cost_seals/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2));
Final_cost_pump(1,1) = C_bop_n_pump(sol(1),sol(2)); % in kW
Final_cost_pcs(1,1) = Cost_pcs;
Final_cost_transformer(1,1) = Cost_transformer;
Final_cost_interconnection(1,1) = Cost_interconnection;
Final_number_of_cells(1,1) = sol(2);
Final_channel_height(1,1) = sol(1);
Final_ASR(1,1) = ASR(sol(1),sol(2));
Final_eta_sys(1,1) = eta_sys(sol(1),sol(2))


clearvars -except Final_cost_sep Final_cost_collector Final_cost_seals Final_cost_pump Final_cost_pcs Final_cost_transformer Final_cost_interconnection Final_number_of_cells Final_channel_height Final_ASR exitflag Final_eta_sys

load ('Catholyte 4p2CB 11MnO2.mat');
load ('Cost components for power_future_mid');
load ('Reciprocating pump parameters');

n_cell_pump = 100;

%% Channel dimensions
channel_width = 2e-2; % 2cm

%% System requirement
power_output = 1000e3; %in W

%% Separator
ASR_sep = 0.1e-4; %R_sep 0.1 ohm.cm2 Hopkins et al.


%% Calculation

% Area specific impedance (in ohm.m2)
ASR_CT_pos = @(channel_height,n) r_ct_pos / (vol_frac_am_pos * channel_height * a_pos);
% ASR_CT_neg = @(channel_height,n) r_ct_neg / (vol_frac_am_neg * channel_height * a_neg);
ASR_CT_neg = 0;
ASR_CONT_pos = @(channel_height,n) r_cont_pos / (vol_frac_am_pos * a_pos);
ASR_CONT_neg = 0;
% ASR_CONT_neg = @(channel_height,n) r_cont_neg / (vol_frac_am_neg * a_neg);

ASR_OHM_pos = @(channel_height,n) (channel_height/2)*(1/catholyte_econd + 1/(catholyte_icond*vol_frac_elec_neg^1.5));
ASR_OHM_neg = @(channel_height,n) (channel_height/2)*(1/anolyte_econd + 1/(anolyte_icond*vol_frac_elec_neg^1.5));

ASR =@(channel_height,n) ASR_sep...
    + ASR_CT_pos(channel_height,n) + ASR_CT_neg...
    + ASR_CONT_pos(channel_height,n) + ASR_CONT_neg...
    + ASR_OHM_pos(channel_height,n) + ASR_OHM_neg(channel_height,n);


% Calculating A_s and delta_pressure_catholyte and anolyte
P_d =@(channel_height,n) (catholyte_phi_eq^2) * (1-eta_v) * eta_v / ASR(channel_height,n);   %Power density required for 1 cell; ASI is an input
A_s =@(channel_height,n) power_output / P_d(channel_height,n); %total area including multiple cells m2
L =@(channel_height,n) (A_s(channel_height,n) / (channel_width * n));
delta_pressure_catholyte =@(channel_height,n) 2 * (catholyte_tau_y/channel_height) * L(channel_height,n); % pressure drop for one cell
% delta_pressure_anolyte =@(vol_frac_am_pos,channel_height,n) 2 * (anolyte_tau_y/channel_height) * L(vol_frac_am_pos,channel_height,n); % pressure drop for one cell

% Calculating flow rate
flow_rate =@(channel_height,n) n_cell_pump * power_output / (n * catholyte_rho_e);

% Calculating efficiency
eps_pump =@(channel_height,n) 2 * (n/n_cell_pump) * delta_pressure_catholyte(channel_height,n) * flow_rate(channel_height,n) / (power_output * pump_eff);
eta_sys =@(channel_height,n) (1-eps_pump(channel_height,n));

% Calculating other costs
C_stack =@(channel_height,n) A_s(channel_height,n)*(Cost_sep + Cost_collector + Cost_seals) ; % in $/W
C_stack_n =@(channel_height,n) (C_stack(channel_height,n) / (power_output*eta_sys(channel_height,n))) * 1000; % in $/kW

pump_power_inHP = @(channel_height,n) delta_pressure_catholyte(channel_height,n) * flow_rate(channel_height,n) * 0.00134102; % in horsepower HP
Cost_pump =@(channel_height,n) (2753 * pump_power_inHP(channel_height,n)^0.57) * 1.59 * Pump_price_factor; % cost per pump
C_bop_n_pump =@(channel_height,n) (2*(n/n_cell_pump)*Cost_pump(channel_height,n) / power_output) * 1000; % in $/kW


% Calculating cost that scales with power
C_power =@(channel_height,n) Cost_BOP_no_pump + Cost_add + C_bop_n_pump(channel_height,n) + C_stack_n(channel_height,n) ; % in kW


lb = [2.5e-4 1000];
ub = [4e-3 14000];

IntCon = 2;


[sol,fval,exitflag,output] = ga(@(x) C_power(x(1),x(2)), 2, [], [], [], [], lb, ub, @(x) nonlinconstr_SSFB_0830_1(x,pump_power_inHP,pump_power_min_inHP), IntCon);

Final_cost_sep(2,1) = (Cost_sep/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2)); % in kW
Final_cost_collector(2,1) = (Cost_collector/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2));
Final_cost_seals(2,1) = (Cost_seals/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2));
Final_cost_pump(2,1) = C_bop_n_pump(sol(1),sol(2)); % in kW
Final_cost_pcs(2,1) = Cost_pcs;
Final_cost_transformer(2,1) = Cost_transformer;
Final_cost_interconnection(2,1) = Cost_interconnection;
Final_number_of_cells(2,1) = sol(2);
Final_channel_height(2,1) = sol(1);
Final_ASR(2,1) = ASR(sol(1),sol(2));
Final_eta_sys(2,1) = eta_sys(sol(1),sol(2));

clearvars -except Final_cost_sep Final_cost_collector Final_cost_seals Final_cost_pump Final_cost_pcs Final_cost_transformer Final_cost_interconnection Final_number_of_cells Final_channel_height Final_ASR exitflag Final_eta_sys

load ('Catholyte 4p2CB 11MnO2.mat');
load ('Cost components for power_future_lower');
load ('Reciprocating pump parameters');

n_cell_pump = 100;

%% Channel dimensions
channel_width = 2e-2; % 2cm

%% System requirement
power_output = 1000e3; %in W

%% Separator
ASR_sep = 0.1e-4; %R_sep 0.1 ohm.cm2 Hopkins et al.


%% Calculation

% Area specific impedance (in ohm.m2)
ASR_CT_pos = @(channel_height,n) r_ct_pos / (vol_frac_am_pos * channel_height * a_pos);
% ASR_CT_neg = @(channel_height,n) r_ct_neg / (vol_frac_am_neg * channel_height * a_neg);
ASR_CT_neg = 0;
ASR_CONT_pos = @(channel_height,n) r_cont_pos / (vol_frac_am_pos * a_pos);
ASR_CONT_neg = 0;
% ASR_CONT_neg = @(channel_height,n) r_cont_neg / (vol_frac_am_neg * a_neg);

ASR_OHM_pos = @(channel_height,n) (channel_height/2)*(1/catholyte_econd + 1/(catholyte_icond*vol_frac_elec_neg^1.5));
ASR_OHM_neg = @(channel_height,n) (channel_height/2)*(1/anolyte_econd + 1/(anolyte_icond*vol_frac_elec_neg^1.5));

ASR =@(channel_height,n) ASR_sep...
    + ASR_CT_pos(channel_height,n) + ASR_CT_neg...
    + ASR_CONT_pos(channel_height,n) + ASR_CONT_neg...
    + ASR_OHM_pos(channel_height,n) + ASR_OHM_neg(channel_height,n);


% Calculating A_s and delta_pressure_catholyte and anolyte
P_d =@(channel_height,n) (catholyte_phi_eq^2) * (1-eta_v) * eta_v / ASR(channel_height,n);   %Power density required for 1 cell; ASI is an input
A_s =@(channel_height,n) power_output / P_d(channel_height,n); %total area including multiple cells m2
L =@(channel_height,n) (A_s(channel_height,n) / (channel_width * n));
delta_pressure_catholyte =@(channel_height,n) 2 * (catholyte_tau_y/channel_height) * L(channel_height,n); % pressure drop for one cell
% delta_pressure_anolyte =@(vol_frac_am_pos,channel_height,n) 2 * (anolyte_tau_y/channel_height) * L(vol_frac_am_pos,channel_height,n); % pressure drop for one cell

% Calculating flow rate
flow_rate =@(channel_height,n) n_cell_pump * power_output / (n * catholyte_rho_e);

% Calculating efficiency
eps_pump =@(channel_height,n) 2 * (n/n_cell_pump) * delta_pressure_catholyte(channel_height,n) * flow_rate(channel_height,n) / (power_output * pump_eff);
eta_sys =@(channel_height,n) (1-eps_pump(channel_height,n));

% Calculating other costs
C_stack =@(channel_height,n) A_s(channel_height,n)*(Cost_sep + Cost_collector + Cost_seals) ; % in $/W
C_stack_n =@(channel_height,n) (C_stack(channel_height,n) / (power_output*eta_sys(channel_height,n))) * 1000; % in $/kW

pump_power_inHP = @(channel_height,n) delta_pressure_catholyte(channel_height,n) * flow_rate(channel_height,n) * 0.00134102; % in horsepower HP
Cost_pump =@(channel_height,n) (2753 * pump_power_inHP(channel_height,n)^0.57) * 1.59 * Pump_price_factor; % cost per pump
C_bop_n_pump =@(channel_height,n) (2*(n/n_cell_pump)*Cost_pump(channel_height,n) / power_output) * 1000; % in $/kW


% Calculating cost that scales with power
C_power =@(channel_height,n) Cost_BOP_no_pump + Cost_add + C_bop_n_pump(channel_height,n) + C_stack_n(channel_height,n) ; % in kW


lb = [2.5e-4 1000];
ub = [4e-3 14000];

IntCon = 2;


[sol,fval,exitflag,output] = ga(@(x) C_power(x(1),x(2)), 2, [], [], [], [], lb, ub, @(x) nonlinconstr_SSFB_0830_1(x,pump_power_inHP,pump_power_min_inHP), IntCon);

Final_cost_sep(3,1) = (Cost_sep/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2)); % in kW
Final_cost_collector(3,1) = (Cost_collector/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2));
Final_cost_seals(3,1) = (Cost_seals/(Cost_sep + Cost_collector + Cost_seals))*C_stack_n(sol(1),sol(2));
Final_cost_pump(3,1) = C_bop_n_pump(sol(1),sol(2)); % in kW
Final_cost_pcs(3,1) = Cost_pcs;
Final_cost_transformer(3,1) = Cost_transformer;
Final_cost_interconnection(3,1) = Cost_interconnection;
Final_number_of_cells(3,1) = sol(2);
Final_channel_height(3,1) = sol(1);
Final_ASR(3,1) = ASR(sol(1),sol(2));
Final_eta_sys(3,1) = eta_sys(sol(1),sol(2));

clearvars -except Final_cost_sep Final_cost_collector Final_cost_seals Final_cost_pump Final_cost_pcs Final_cost_transformer Final_cost_interconnection Final_number_of_cells Final_channel_height Final_ASR exitflag Final_eta_sys

T = table(Final_channel_height,Final_number_of_cells,Final_ASR,Final_cost_sep,Final_cost_collector,Final_cost_seals,Final_cost_pump,Final_cost_pcs,Final_cost_transformer,Final_cost_interconnection,Final_eta_sys)
