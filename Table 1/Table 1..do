**Column 1 Sequeira 2016
reg lba  tariff_change_post2008 tariff_change_2008  tariff2007 differentiated agri perishable dfs i.clear_agent  lvalue_tonnage day_w_arrival  psi monitor post_2008  i.hc_group hc_4digits rsa term  diff_post_2008 agri_post_2008 lvalue_ton_post_2008  perishable_post_2008 day_w_arrival_post2008 dfs_post_2008 psi_post_2008 tariff2007_post_2008, vce(cluster hc_4digits)


**Column 2
gen inter_diff=tariff_change_post2008*differentiated
gen inter_agri=tariff_change_post2008*agri
gen inter_psi=tariff_change_post2008*psi
gen inter_perish=tariff_change_post2008*perishable
gen inter_dfs=tariff_change_post2008*dfs
gen inter_lvalue=tariff_change_post2008*lvalue_tonnage
gen inter_dwa=tariff_change_post2008*day_w_arrival
gen inter_rsa=tariff_change_post2008*rsa
gen inter_tariff=tariff_change_post2008*tariff2007
gen inter_monitor=tariff_change_post2008*monitor



reg lba  tariff_change_post2008 inter_diff inter_agri inter_lvalue inter_perish inter_dwa inter_dfs inter_psi inter_rsa inter_tariff inter_monitor tariff_change_2008  tariff2007 differentiated agri perishable dfs i.clear_agent  lvalue_tonnage day_w_arrival  psi monitor post_2008  i.hc_group hc_4digits rsa term  diff_post_2008 agri_post_2008 lvalue_ton_post_2008  perishable_post_2008 day_w_arrival_post2008 dfs_post_2008 psi_post_2008 tariff2007_post_2008, vce(cluster hc_4digits)
