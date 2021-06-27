
% Derived constants

molarity_am_pos = vol_frac_am_pos * (rho_pos_am) / molecular_weight_pos_am; %in mol/L
catholyte_rho_e = no_of_elec_pos * molarity_am_pos * 1000 * F * catholyte_phi_eq * eta_v; %in J/m3
anolyte_rho_e = catholyte_rho_e; %in J/m3
molarity_am_neg = anolyte_rho_e/(no_of_elec_neg * F * 1000 * catholyte_phi_eq * eta_v); %in mol/L
vol_frac_am_neg = molarity_am_neg * molecular_weight_neg_am / (rho_neg_am); % <1

vol_frac_elec_pos = 1-vol_frac_poly_pos-vol_frac_condadd_pos-vol_frac_am_pos;
vol_frac_elec_neg = 1-vol_frac_poly_neg-vol_frac_condadd_neg-vol_frac_am_neg;


x_neg = [vol_frac_elec_neg; vol_frac_poly_neg; vol_frac_condadd_neg; vol_frac_am_neg];
rho_neg = [rho_neg_elec; rho_neg_poly; rho_neg_condadd; rho_neg_am];
c_neg = [c_neg_elec; c_neg_poly; c_neg_condadd; c_neg_am];

x_pos = [vol_frac_elec_pos; vol_frac_poly_pos; vol_frac_condadd_pos; vol_frac_am_pos];
rho_pos = [rho_pos_elec; rho_pos_poly; rho_pos_condadd; rho_pos_am];
c_pos = [c_pos_elec; c_pos_poly; c_pos_condadd; c_pos_am];

C_am_neg = sum(x_neg .* rho_neg .* c_neg);
C_am_pos = sum(x_pos .* rho_pos .* c_pos);

% Total am cost (in $/kWh)
C_am = (C_am_pos/(catholyte_rho_e/3600) + C_am_neg/(anolyte_rho_e/3600))*1000;

% Total overnight energy cost ($/kWh)
C_tank = (c_tank/(catholyte_rho_e/3600) + c_tank/(anolyte_rho_e/3600))*1000;
