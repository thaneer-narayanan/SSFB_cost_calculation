function Capital_cost

clear all

%%
clear all

load('Catholyte 4p2CB 11MnO2_upper')
Final_CC_energy(1,1) = C_energy;
load('Catholyte 4p2CB 11MnO2')
Final_CC_energy(2,1) = C_energy;
load('Catholyte 4p2CB 11MnO2_lower')
Final_CC_energy(3,1) = C_energy;

clearvars -except Final_CC_energy 

load('!Results_future_4p2CB 11MnO2_n_cell_pump_100')

Final_CC_stack = Final_cost_collector + Final_cost_sep + Final_cost_seals;
Final_CC_BOP = Final_cost_pump + Final_cost_pcs + Final_cost_hex + Final_cost_transformer + Final_cost_interconnection;

Final_CC_power = Final_CC_stack + Final_CC_BOP;
Final_CC_energy = Final_CC_energy ./ Final_eta_sys;

Cost_add = 197;
f_install = 0.2;
duration = logspace(0,3,30);

Final_CC = zeros(length(Final_CC_power),length(duration));
for i=1:length(Final_CC_power)
    for j=1:length(duration)
        Final_CC(i,j) = (Final_CC_power(i) / duration(j) + Final_CC_energy(i)) * (1+f_install) + Cost_add / duration(j);
    end
end

clearvars -except Final_CC_energy Final_CC_stack Final_CC_BOP Final_CC_power Cost_add f_install duration Final_CC

duration = transpose(duration);
Final_CC = transpose(Final_CC);
T = table(duration,Final_CC(:,1),Final_CC(:,2),Final_CC(:,3));

filename = '!Results_capitalcost_future_4p2CB 11MnO2_n_cell_pump_100.mat';
save(filename)

