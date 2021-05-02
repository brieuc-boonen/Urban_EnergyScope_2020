# -------------------------------------------------------------------------------------------------------------------------													
#	EnergyScope TD is an open-source energy model suitable for country scale analysis. It is a simplified representation of an urban or national energy system accounting for the energy flows												
#	within its boundaries. Based on a hourly resolution, it optimises the design and operation of the energy system while minimizing the cost of the system.												
#													
#	Copyright (C) <2018-2019> <Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland and Université catholique de Louvain (UCLouvain), Belgium>
#													
#	Licensed under the Apache License, Version 2.0 (the "License");												
#	you may not use this file except in compliance with the License.												
#	You may obtain a copy of the License at												
#													
#		http://www.apache.org/licenses/LICENSE-2.0												
#													
#	Unless required by applicable law or agreed to in writing, software												
#	distributed under the License is distributed on an "AS IS" BASIS,												
#	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.												
#	See the License for the specific language governing permissions and												
#	limitations under the License.												
#													
#	Description and complete License: see LICENSE file.												
# -------------------------------------------------------------------------------------------------------------------------		


######################### 
###  SETS [Figure 3]  ###
#########################

## MAIN SETS: Sets whose elements are input directly in the data file
set PERIODS := 1 .. 8760; # time periods (hours of the year)
set HOURS := 1 .. 24; # hours of the day
set TYPICAL_DAYS:= 1 .. 12; # typical days
set T_H_TD within {PERIODS, HOURS, TYPICAL_DAYS}; # set linking periods, hours, days, typical days
set SECTORS; # sectors of the energy system
set END_USES_INPUT; # Types of demand (end-uses). Input to the model
set END_USES_CATEGORIES; # Categories of demand (end-uses): electricity, heat, mobility
set END_USES_TYPES_OF_CATEGORY {END_USES_CATEGORIES}; # Types of demand (end-uses).
set RESOURCES; # Resources: fuels (renewables and fossils) and electricity imports
set BIOFUELS within RESOURCES; # imported biofuels.
set EXPORT within RESOURCES; # exported resources
set IMPORT_FEED_IN within RESOURCES; 
set ELECTRANSFER within RESOURCES;
set END_USES_TYPES := setof {i in END_USES_CATEGORIES, j in END_USES_TYPES_OF_CATEGORY [i]} j; # secondary set
set TECHNOLOGIES_OF_END_USES_TYPE {END_USES_TYPES}; # set all energy conversion technologies (excluding storage technologies and infrastructure)
set HOME_TECHNOLOGIES {END_USES_TYPES};
set STORAGE_TECH; #  set of storage technologies 
set STORAGE_OF_END_USES_TYPES {END_USES_TYPES} within STORAGE_TECH; # set all storage technologies related to an end-use types (used for thermal solar (TS))
set INFRASTRUCTURE; # Infrastructure: DHN, grid, and intermediate energy conversion technologies (i.e. not directly supplying end-use demand)
set RENOVATION;
set REFSIZETECH;

## SECONDARY SETS: a secondary set is defined by operations on MAIN SETS
set LAYERS := (RESOURCES diff BIOFUELS diff EXPORT diff IMPORT_FEED_IN) union END_USES_TYPES; # Layers are used to balance resources/products in the system
set TECHNOLOGIES := (setof {i in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE [i]} j) union STORAGE_TECH union INFRASTRUCTURE union RENOVATION; 
set TECHNOLOGIES_OF_END_USES_CATEGORY {i in END_USES_CATEGORIES} within TECHNOLOGIES := setof {j in END_USES_TYPES_OF_CATEGORY[i], k in TECHNOLOGIES_OF_END_USES_TYPE [j]} k;
set RE_RESOURCES within RESOURCES; # List of RE resources (including wind hydro solar), used to compute the RE share
set V2G within TECHNOLOGIES;   # EVs which can be used for vehicle-to-grid (V2G).
set EVs_BATT   within STORAGE_TECH; # specific battery of EVs
set EVs_BATT_OF_V2G {V2G}; # Makes the link between batteries of EVs and the V2G technology
set STORAGE_DAILY within STORAGE_TECH;# Storages technologies for daily application 
set TS_OF_DEC_TECH {TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"} diff {"DEC_HP_ELEC_A2A"} diff {"DEC_STOVE_WOOD"}} ; # Makes the link between TS and the technology producing the heat

##Additional SETS added just to simplify equations.
set TYPICAL_DAY_OF_PERIOD {t in PERIODS} := setof {h in HOURS, td in TYPICAL_DAYS: (t,h,td) in T_H_TD} td; #TD_OF_PERIOD(T)
set HOUR_OF_PERIOD {t in PERIODS} := setof {h in HOURS, td in TYPICAL_DAYS: (t,h,td) in T_H_TD} h; #H_OF_PERIOD(T)

## Additional SETS: only needed for printing out results (not represented in Figure 3).
set COGEN within TECHNOLOGIES; # cogeneration tech;
set BOILERS within TECHNOLOGIES; # boiler tech
set HPS within TECHNOLOGIES;

#################################
### PARAMETERS [Tables 1-2]   ###
#################################

## Parameters added to include time series in the model [Table 1]
param electricity_time_series {HOURS, TYPICAL_DAYS} >= 0, <= 1; # %_elec [-]: factor for sharing lighting across typical days (adding up to 1)
param heating_time_series {HOURS, TYPICAL_DAYS} >= 0, <= 1; # %_sh [-]: factor for sharing space heating across typical days (adding up to 1)
param mob_pass_time_series {HOURS, TYPICAL_DAYS} >= 0, <= 1; # %_pass [-]: factor for sharing passenger transportation across Typical days (adding up to 1) based on https://www.fhwa.dot.gov/policy/2013cpr/chap1.cfm
param mob_freight_time_series {HOURS, TYPICAL_DAYS} >= 0, <= 1; # %_fr [-]: factor for sharing freight transportation across Typical days (adding up to 1)
param c_p_t {TECHNOLOGIES, HOURS, TYPICAL_DAYS} default 1; #Hourly capacity factor [-]. If = 1 (default value) <=> no impact.
param cop_time_series {RESOURCES union TECHNOLOGIES union STORAGE_TECH, LAYERS, HOURS, TYPICAL_DAYS} default 1;

## Parameters added to define scenarios and technologies [Table 2]
param end_uses_demand_year {END_USES_INPUT, SECTORS} >= 0 default 0; # end_uses_year [GWh]: table end-uses demand vs sectors (input to the model). Yearly values. [Mpkm] or [Mtkm] for passenger or freight mobility.
param end_uses_input {i in END_USES_INPUT} := sum {s in SECTORS} (end_uses_demand_year [i,s]); # end_uses_input (Figure 1.4) [GWh]: total demand for each type of end-uses across sectors (yearly energy) as input from the demand-side model. [Mpkm] or [Mtkm] for passenger or freight mobility.
param auto_consumption_rate >= 0;
param auto_sufficiancy_rate >= -10;
param i_rate > 0; # discount rate [-]: real discount rate
param re_share_primary >= 0; # re_share [-]: minimum share of primary energy coming from RE*/

#param TotalGWP >= 0; # GWP_tot [ktCO2-eq./year]: Total global warming potential (GWP) emissions in the system
param gwp_limit >= 0;    # [ktCO2-eq./year] maximum gwp emissions allowed.
param share_mobility_public_min >= 0, <= 1; # %_public,min [-]: min limit for penetration of public mobility over total mobility 
param share_mobility_public_max >= 0, <= 1; # %_public,max [-]: max limit for penetration of public mobility over total mobility 
param share_renovation_min >= 0, <= 1; # %_public,min [-]: min limit for penetration of public mobility over total mobility 
param share_renovation_max >= 0, <= 1; # %_public,max [-]: max limit for penetration of public mobility over total mobility 
param share_freight_train_min >= 0, <= 1; # %_rail,min [-]: min limit for penetration of train in freight transportation
param share_freight_train_max >= 0, <= 1; # %_rail,min [-]: max limit for penetration of train in freight transportation
param share_heat_dhn_min >= 0, <= 1; # %_dhn,min [-]: min limit for penetration of dhn in low-T heating
param share_heat_dhn_max >= 0, <= 1; # %_dhn,max [-]: max limit for penetration of dhn in low-T heating
param t_op {HOURS, TYPICAL_DAYS} default 1;# [h]: operating time 
param f_max {TECHNOLOGIES} >= 0; # Maximum feasible installed capacity [GW], refers to main output. storage level [GWh] for STORAGE_TECH
param f_min {TECHNOLOGIES} >= 0; # Minimum feasible installed capacity [GW], refers to main output. storage level [GWh] for STORAGE_TECH
param fmax_perc {TECHNOLOGIES} >= 0, <= 1 default 1; # value in [0,1]: this is to fix that a technology can at max produce a certain % of the total output of its sector over the entire year
param fmin_perc {TECHNOLOGIES} >= 0, <= 1 default 0; # value in [0,1]: this is to fix that a technology can at min produce a certain % of the total output of its sector over the entire year
param avail {RESOURCES} >= 0; # Yearly availability of resources [GWh/y]
param c_op {RESOURCES} >= -1; # cost of resources in the different periods [MCHF/GWh]
param n_car_max >=0; #  [car] Maximum amount of cars. Required to compute the aggregated size of EVs batteries.
param peak_sh_factor >= 0;   # %_Peak_sh [-]: ratio between highest yearly demand and highest TDs demand
param peak_real_factor >= 0;   # 
param worst_year_factor >= 0;   # 
param layers_in_out {RESOURCES union TECHNOLOGIES diff STORAGE_TECH , LAYERS}; # f: input/output Resources/Technologies to Layers. Reference is one unit ([GW] or [Mpkm/h] or [Mtkm/h]) of (main) output of the resource/technology. input to layer (output of technology) > 0.
param ref_size {TECHNOLOGIES} >= 0; # f_ref: reference size of each technology, expressed in the same units as the layers_in_out table. Refers to main output (heat for cogen technologies). storage level [GWh] for STORAGE_TECH
param c_inv {TECHNOLOGIES} >= 0; # Specific investment cost [MCHF/GW].[MCHF/GWh] for STORAGE_TECH
param c_maint {TECHNOLOGIES} >= 0; # O&M cost [MCHF/GW/year]: O&M cost does not include resource (fuel) cost. [MCHF/GWh/year] for STORAGE_TECH
param lifetime {TECHNOLOGIES} >= 0; # n: lifetime [years]
param tau {i in TECHNOLOGIES} := i_rate * (1 + i_rate)^lifetime [i] / (((1 + i_rate)^lifetime [i]) - 1); # Annualisation factor ([-]) for each different technology [Eq. 2]
param gwp_constr {TECHNOLOGIES} >= 0; # GWP emissions associated to the construction of technologies [ktCO2-eq./GW]. Refers to [GW] of main output
param gwp_op {RESOURCES} >= 0; # GWP emissions associated to the use of resources [ktCO2-eq./GWh]. Includes extraction/production/transportation and combustion
param c_p {TECHNOLOGIES} >= 0, <= 1 default 1; # yearly capacity factor of each technology [-], defined on annual basis. Different than 1 if sum {t in PERIODS} F_t (t) <= c_p * F
param storage_eff_in {STORAGE_TECH , LAYERS} >= 0, <= 1; # eta_sto_in [-]: efficiency of input to storage from layers.  If 0 storage_tech/layer are incompatible
param storage_eff_out {STORAGE_TECH , LAYERS} >= 0, <= 1; # eta_sto_out [-]: efficiency of output from storage to layers. If 0 storage_tech/layer are incompatible
param storage_losses {STORAGE_TECH} >= 0, <= 1; # %_sto_loss [-]: Self losses in storage (required for Li-ion batteries). Value = self discharge in 1 hour.
param storage_charge_time    {STORAGE_TECH} >= 0; # t_sto_in [h]: Time to charge storage (Energy to Power ratio). If value =  5 <=>  5h for a full charge.
param storage_discharge_time {STORAGE_TECH} >= 0; # t_sto_out [h]: Time to discharge storage (Energy to Power ratio). If value =  5 <=>  5h for a full discharge.
param storage_availability {STORAGE_TECH} >=0, default 1;# %_sto_avail [-]: Storage technology availability to charge/discharge. Used for EVs 
param loss_network {END_USES_TYPES} >= 0 default 0; # %_net_loss: Losses coefficient [0; 1] in the networks (grid and DHN)
param Batt_per_Car {V2G} >= 0; # ev_Batt_size [GWh]: Battery size per EVs car technology
param c_grid_extra >=0; # Cost to reinforce the grid due to IRE penetration [MCHF].
param mult_factor {TECHNOLOGIES} >=0 default 0; # Multiplicative factor that is applied to the tradable green certificates (TGC) to give more (>1) or less (<1) value to the traded GC for each tech.
param c_gc >=0; # Mean price of one TGC in today's market
param living_area >0; # Living area in [m2]
param c_inc_fix {TECHNOLOGIES} >=0 default 0; # Incentive revenue
param c_inc_var {TECHNOLOGIES} >=0 default 0; # Incentive revenue
##Additional parameter (not presented in the paper)
param total_time := sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} (t_op [h, td]); # [h]. added just to simplify equations
param roof_limit default 6; # 
param ngroup default 3; #
param pop default 1; #
param lin_dens default 2;
param tax default 1;
param carbon_tax default 0;


#################################
###   VARIABLES [Tables 3-4]  ###
#################################


##Independent variables [Table 3] :
var Share_Mobility_Public >= share_mobility_public_min, <= share_mobility_public_max; # %_Public: Ratio [0; 1] public mobility over total passenger mobility
var Share_Renovation >= share_renovation_min, <= share_renovation_max;
var Share_Freight_Train, >= share_freight_train_min, <= share_freight_train_max; # %_Rail: Ratio [0; 1] rail transport over total freight transport
var Share_Heat_Dhn, >= share_heat_dhn_min, <= share_heat_dhn_max; # %_DHN: Ratio [0; 1] centralized over total low-temperature heat
var Share_Heat_Dec, >= 0, <= 1;
var Share_Heat_Ren, >= 0, <= 1;
var F {TECHNOLOGIES} >= 0; # F: Installed capacity ([GW]) with respect to main output (see layers_in_out). [GWh] for STORAGE_TECH.
var F_t {RESOURCES union TECHNOLOGIES, HOURS, TYPICAL_DAYS} >= 0; # F_t: Operation in each period [GW] or, for STORAGE_TECH, storage level [GWh]. multiplication factor with respect to the values in layers_in_out table. Takes into account c_p
var Storage_in {i in STORAGE_TECH, LAYERS, HOURS, TYPICAL_DAYS} >= 0; # Sto_in [GW]: Power input to the storage in a certain period
var Storage_out {i in STORAGE_TECH, LAYERS, HOURS, TYPICAL_DAYS} >= 0; # Sto_out [GW]: Power output from the storage in a certain period
var Shares_Mobility_Passenger {TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_PASSENGER"]} >=0; # %_MobPass [-]: Constant share of passenger mobility
var Shares_LowT_Dec {TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}}>=0 ; # %_HeatDec [-]: Constant share of heat Low T decentralised + its specific thermal solar
var Shares_LowT_DHN >=0 ; # %_HeatDec [-]: Constant share of heat Low T decentralised + its specific thermal solar
var F_Solar         {TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] union TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DEC_SOLAR"} diff {"DHN_SOLAR"}} >=0; # F_sol [GW]: Solar thermal installed capacity per heat decentralised technologies
var F_t_Solar       {TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] union TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DEC_SOLAR"} diff {"DHN_SOLAR"}, h in HOURS, td in TYPICAL_DAYS} >= 0; # F_t_sol [GW]: Solar thermal operating per heat decentralised technologies
var X_active {TECHNOLOGIES} binary default 0;

##Dependent variables [Table 4] :
var End_Uses {LAYERS, HOURS, TYPICAL_DAYS} >= 0; #EndUses [GW]: total demand for each type of end-uses (hourly power). Defined for all layers (0 if not demand). [Mpkm] or [Mtkm] for passenger or freight mobility.
var Number_Of_Units {TECHNOLOGIES} integer; # N: number of units of size ref_size which are installed.
var TotalCost >= -100000; # C_tot [ktCO2-eq./year]: Total GWP emissions in the system.
var C_inv {TECHNOLOGIES} >= 0; #C_inv [MCHF]: Total investment cost of each technology
var C_maint {TECHNOLOGIES} >= 0; #C_maint [MCHF/year]: Total O&M cost of each technology (excluding resource cost)
var C_op {RESOURCES} >= -100000; #C_op [MCHF/year]: Total O&M cost of each resource
var TotalGWP >= 0; # GWP_tot [ktCO2-eq./year]: Total global warming potential (GWP) emissions in the system
var GWP_constr {TECHNOLOGIES} >= 0; # GWP_constr [ktCO2-eq.]: Total emissions of the technologies
var GWP_op {RESOURCES} >= 0; #  GWP_op [ktCO2-eq.]: Total yearly emissions of the resources [ktCO2-eq./y]
var Network_losses {END_USES_TYPES, HOURS, TYPICAL_DAYS} >= 0; # Net_loss [GW]: Losses in the networks (normally electricity grid and DHN)
var Storage_level {STORAGE_TECH, PERIODS} >= 0; # Sto_level [GWh]: Energy stored at each period
var Prosumer_tax >=0;
#var R_gc  {TECHNOLOGIES} >=0 default 0;
var R_inc {TECHNOLOGIES} >=0 default 0;
var R_ren {TECHNOLOGIES} >=0 default 0;

#param GWP_op {RESOURCES} = 0;

#########################################
###      CONSTRAINTS Eqs [1-42]       ###
#########################################

## End-uses demand calculation constraints 
#-----------------------------------------

# [Figure 4] From annual energy demand to hourly power demand. End_Uses is non-zero only for demand layers.
subject to end_uses_t {l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	End_Uses [l, h, td] = (if l == "ELECTRICITY" 
		then
			(end_uses_input[l] / total_time + end_uses_input["LIGHTING"] * electricity_time_series [h, td] / t_op [h, td] ) * pop + Network_losses [l,h,td]
		else (if l == "HEAT_LOW_T_DHN" then
			(pop * end_uses_input["HEAT_LOW_T_HW"] / total_time + (pop * end_uses_input ["HEAT_LOW_T_SH"] ) * heating_time_series [h, td] / t_op [h, td] ) * Share_Heat_Dhn + Network_losses [l,h,td]
		else (if l == "HEAT_LOW_T_DECEN" then
			(pop * end_uses_input["HEAT_LOW_T_HW"] / total_time + (pop * end_uses_input ["HEAT_LOW_T_SH"] ) * heating_time_series [h, td] / t_op [h, td] ) * Share_Heat_Dec
		else (if l == "MOB_PUBLIC" then
			(end_uses_input["MOBILITY_PASSENGER"] * mob_pass_time_series [h, td] / t_op [h, td]  ) * pop * Share_Mobility_Public
		else (if l == "MOB_PRIVATE" then
			(end_uses_input["MOBILITY_PASSENGER"] * mob_pass_time_series [h, td] / t_op [h, td]  ) * pop * (1 - Share_Mobility_Public)
		else (if l == "MOB_FREIGHT_RAIL" then
			(end_uses_input["MOBILITY_FREIGHT"]   * mob_freight_time_series [h, td] / t_op [h, td] ) * pop *  Share_Freight_Train
		else (if l == "MOB_FREIGHT_ROAD" then
			(end_uses_input["MOBILITY_FREIGHT"]   * mob_freight_time_series [h, td] / t_op [h, td] ) * pop * (1 - Share_Freight_Train)
		else 
			0 ))))))); # For all layers which don't have an end-use demand


subject to Share_RENOVATION_total : 
Share_Renovation = (sum{k in RENOVATION} (F [k] / ( pop*(end_uses_input ["HEAT_LOW_T_SH"] + end_uses_input["HEAT_LOW_T_HW"] ))));

subject to ShareDHN_implementation_LowTH : # We force these terms to equal 1 <=> satisfy the total demand
	(Share_Heat_Dhn + Share_Heat_Dec + Share_Heat_Ren) = 1;

subject to renovation_share :
	Share_Heat_Ren = sum{k in RENOVATION} (F [k] / ( pop*(end_uses_input ["HEAT_LOW_T_SH"] + end_uses_input["HEAT_LOW_T_HW"] ))) ;

# Multiplication factor

# [Eq. 1.7] Number of purchased technologies. Integer variable (so that we have only integer multiples of the reference size)
/*subject to number_of_units {j in REFSIZETECH}:
	Number_Of_Units [j] = F [j] / ref_size [j];*/

## Cost
#------

# [Eq. 1]
subject to totalcost_cal:
TotalCost = (sum {j in TECHNOLOGIES diff RENOVATION} (tau [j] * ( C_inv [j] /*- R_gc [j]*/ - R_inc [j] ) + C_maint [j] ) + sum {k in RENOVATION} ( tau [k] * (C_inv [k] - R_ren [k]) ) + sum {i in RESOURCES} C_op [i]) * tax /*+ Prosumer_tax*/; 

# [Eq. 3] Investment cost of each technology
subject to investment_cost_calc {j in TECHNOLOGIES}: 
	C_inv [j] = c_inv [j] * F [j];
		
# [Eq. 4] O&M cost of each technology
subject to main_cost_calc {j in TECHNOLOGIES}: 
	C_maint [j] = c_maint [j] * F [j];		

# [Eq. 5] Total cost of each resource
subject to op_cost_calc {i in RESOURCES}:
	C_op [i] = sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ((c_op [i] + gwp_op [i]*carbon_tax) * F_t [i, h, td] * t_op [h, td]) ;

# [PoC] Total revenue with Green Certificates
/*
subject to rgc_cost_calc {j in TECHNOLOGIES}:
	R_gc [j] = c_gc * mult_factor[j] * sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ( F_t [j, h, td] * t_op [h, td] );
*/

# [PoC] Total revenue with other incentives
subject to rinc_cost_cond {j in TECHNOLOGIES diff RENOVATION}:
	(0.7 * C_inv [j]) >= R_inc [j] ;  /*+ c_inc_var [j] * F[j]*/ #sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ( F_t [j, h, td] * t_op [h, td] );
subject to rinc_cost_calc {j in TECHNOLOGIES diff RENOVATION}:
	R_inc [j] <= (c_inc_fix [j] * ngroup * 0 /* X_active [j]  /*+ c_inc_var [j] * F[j]*/); #sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ( F_t [j, h, td] * t_op [h, td] );

subject to rren_cost_cond {j in RENOVATION}:
	(0.7 * C_inv [j]) >= R_ren [j] ; 
subject to rren_cost_calc {j in RENOVATION}:
	R_ren [j] <= (c_inc_var [j] * ngroup * F[j]);
	
## Emissions
#-----------

# [Eq. 6]
subject to totalGWP_calc:
	TotalGWP = sum {j in TECHNOLOGIES} (GWP_constr [j] / lifetime [j]) + sum {i in RESOURCES} GWP_op [i];
	
# [Eq. 7]
subject to gwp_constr_calc {j in TECHNOLOGIES}:
	GWP_constr [j] = gwp_constr [j] * F [j];

# [Eq. 8]
subject to gwp_op_calc {i in RESOURCES}:
	GWP_op [i] = gwp_op [i] * sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ( F_t [i, h, td] * t_op [h, td] );	


## Multiplication factor
#-----------------------

# [Eq. 9] min & max limit to the size of each technology

/*
subject to size_limit_min {j in TECHNOLOGIES diff TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"] diff INFRASTRUCTURE}:
	F [j] >= f_min [j] * pop * X_active [j];
subject to size_limit_max {j in TECHNOLOGIES diff TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"] diff INFRASTRUCTURE}:
	F [j] <= f_max [j] * pop * X_active [j];
*/
/*
subject to size_limit_min_dhn {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] union STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"] union INFRASTRUCTURE}:
	F [j] >= f_min [j] * X_active [j];
*/





/*
subject to renov_install :
	F ["FACADE"] = 12937680 ;
*/

subject to size_limit_max_rest {j in TECHNOLOGIES diff TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"]}:
	F [j] <= f_max [j] * pop /*X_active [j]*/;

/*subject to size_limit_min_rest {j in TECHNOLOGIES diff TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"] diff INFRASTRUCTURE}:
	F [j] >= f_min [j] * pop * X_active [j];*/

subject to size_limit_max_rest_dhn {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] union STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"] union TECHNOLOGIES_OF_END_USES_TYPE["ELECTRICITY"]}:
	F [j] <= f_max [j] /*X_active [j]*/;
	
/*subject to size_limit_min_dhn_dhn {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] union STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"] union INFRASTRUCTURE}:
	F [j] >= f_min [j] * X_active [j];*/

/*subject to size_limit_min_hp:
	F["DHN_BOILER_WOOD_S"] >= f_min["DHN_BOILER_WOOD_S"];*/
/*subject to size_limit_act_rest {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] union STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DECEN"]}:
	X_active [j] = 0 ;*/

	
# [Eq. 10] relation between power and capacity via period capacity factor. This forces max hourly output (e.g. renewables)
subject to capacity_factor_t {j in TECHNOLOGIES, h in HOURS, td in TYPICAL_DAYS}:
	F_t [j, h, td] <= F [j] * c_p_t [j, h, td] ;
	
# [Eq. 11] relation between mult_t and mult via yearly capacity factor. This one forces total annual output
subject to capacity_factor {j in TECHNOLOGIES}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} (F_t [j, h, td] * t_op [h, td]) <= F [j] * c_p [j] * total_time;

## Resources
#-----------

# [Eq. 12] Resources availability equation
subject to resource_availability {i in RESOURCES}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [i, h, td] * t_op [h, td]) <= pop * avail [i];

## Layers
#--------

# [Eq. 13] Layer balance equation with storage. Layers: input > 0, output < 0. Demand > 0. Storage: in > 0, out > 0;
# output from technologies/resources/storage - input to technologies/storage = demand. Demand has default value of 0 for layers which are not end_uses
subject to layer_balance {l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
		sum {i in RESOURCES union TECHNOLOGIES diff STORAGE_TECH } 
		(layers_in_out[i, l] / cop_time_series[i,l, h, td]  * F_t [i, h, td])	
		+ sum {j in STORAGE_TECH} ( Storage_out [j, l, h, td] - Storage_in [j, l, h, td] )
		- End_Uses [l, h, td]
		= 0;

## Storage	
#---------
	
# [Eq. 14] The level of the storage represents the amount of energy stored at a certain time.
subject to storage_level {j in STORAGE_TECH, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]}:
	Storage_level [j, t] = (if t == 1 then
	 			Storage_level [j, card(PERIODS)] * (1.0 -  storage_losses[j])
				+ t_op [h, td] * (   (sum {l in LAYERS: storage_eff_in [j,l] > 0}  (Storage_in [j, l, h, td]  * storage_eff_in  [j, l])) 
				                   - (sum {l in LAYERS: storage_eff_out [j,l] > 0} (Storage_out [j, l, h, td] / storage_eff_out [j, l])))
	else
	 			Storage_level [j, t-1] * (1.0 -  storage_losses[j])
				+ t_op [h, td] * (   (sum {l in LAYERS: storage_eff_in [j,l] > 0}  (Storage_in [j, l, h, td]  * storage_eff_in  [j, l])) 
				                   - (sum {l in LAYERS: storage_eff_out [j,l] > 0} (Storage_out [j, l, h, td] / storage_eff_out [j, l])))
				);

# [Eq. 15] Bounding daily storage
subject to impose_daily_storage {j in STORAGE_DAILY, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]}:
	Storage_level [j, t] = F_t [j, h, td];
	
# [Eq. 16] Bounding seasonal storage
subject to limit_energy_stored_to_maximum {j in STORAGE_TECH diff STORAGE_DAILY , t in PERIODS}:
	Storage_level [j, t] <= F [j];# Never exceed the size of the storage unit
	
# [Eqs. 17-18] Each storage technology can have input/output only to certain layers. If incompatible then the variable is set to 0
subject to storage_layer_in {j in STORAGE_TECH, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_in [j, l, h, td] * (ceil (storage_eff_in [j, l]) - 1) = 0;
subject to storage_layer_out {j in STORAGE_TECH, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_out [j, l, h, td] * (ceil (storage_eff_out [j, l]) - 1) = 0;
		
# [Eq. 19] limit the Energy to power ratio. 
subject to limit_energy_to_power_ratio {j in STORAGE_TECH , l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_in [j, l, h, td] * storage_charge_time[j] + Storage_out [j, l, h, td] * storage_discharge_time[j] <=  F [j] * storage_availability[j];


## Infrastructure
#----------------

# [Eq. 20] Calculation of losses for each end-use demand type (normally for electricity and DHN)
subject to network_losses {eut in END_USES_TYPES, h in HOURS, td in TYPICAL_DAYS}:
	Network_losses [eut,h,td] = (sum {j in RESOURCES union TECHNOLOGIES diff STORAGE_TECH: layers_in_out [j, eut] > 0} ((layers_in_out[j, eut]) * F_t [j, h, td])) * loss_network [eut];

# [Eq. 21] 9.4 BCHF is the extra investment needed if there is a big deployment of stochastic renewables
/*subject to extra_grid:
	F ["GRID"] = 1 + (c_grid_extra / c_inv["GRID"]) * (F ["WIND"] + F ["PV"]) / (f_max ["WIND"] + f_max ["PV"]);*/

/*
subject to lines_dhn :
	F ["DHN_LINES"] = pop * Share_Heat_Dhn;

let lin_dens :=5;
let c_inv["DHN_LINES"] := 1386*(18/lin_dens) + 2175; # UK prices
*/


subject to singlelines_dhn :
	F ["DHN_SINGLE"] =  ( /*(1 - Share_Heat_Ren - Share_Heat_Dec) * */(end_uses_input["HEAT_LOW_T_HW"] + end_uses_input["HEAT_LOW_T_SH"]) ) / 1000 / lin_dens * pop;

subject to servicelines_dhn :
	F ["DHN_SERVICE"] = pop;


# [Eq. 22] DHN: assigning a cost to the network
subject to extra_dhn {k in {"DHN_SUBSTA"} union {"DHN_STATIONS"}}:
	F [k] = Share_Heat_Dhn * max {h in HOURS, td in TYPICAL_DAYS} ((pop * end_uses_input["HEAT_LOW_T_HW"] / total_time + (pop * end_uses_input ["HEAT_LOW_T_SH"] ) * heating_time_series [h, td] / t_op [h, td] ) ); # * Share_Heat_Dhn
	
/*	
subject to price_lines :
	c_inv["DHN_LINES"] = 1386*(18/lin_dens) + 2175 # UK prices
#	c_inv["DHN_LINES"] = 1000*(18/lin_dens) + 3200 # DNK prices
;
*/

	
	#sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}} (2 * F [j]);
	#End_Uses ["HEAT_LOW_T_DHN", h, td];
	#max {h in HOURS, td in TYPICAL_DAYS} ((pop * end_uses_input["HEAT_LOW_T_HW"] / total_time + (pop * end_uses_input ["HEAT_LOW_T_SH"] ) * heating_time_series [h, td] / t_op [h, td] ) * Share_Heat_Dhn);
	#sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}, i in STORAGE_OF_END_USES_TYPES["HEAT_LOW_T_DHN"]} (F [j] + 0.01*F[i]/storage_discharge_time[i]);
#	F ["DHN"] = sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}} (F [j] + F_Solar[j]);

#	F ["DHN"] = sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"], h in HOURS, td in TYPICAL_DAYS} (F_t [j, h, td] * layers_in_out[j, "HEAT_LOW_T_DHN"] / cop_time_series [j,"HEAT_LOW_T_DHN", h, td]);


#subject to thermal_solar_capacity_factor_dhn {h in HOURS, td in TYPICAL_DAYS}:
#	F_t["DHN_SOLAR", h, td] <= sum{j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}} (F[j] * c_p_t["DHN_SOLAR", h, td]);
	

## Additional constraints
#------------------------

/*
# [Eq. 25] Operating strategy in mobility passenger (to make model more realistic)
# Each passenger mobility technology (j) has to supply a constant share  (Shares_Mobility_Passenger[j]) of the passenger mobility demand
subject to operating_strategy_mob_passenger{j in TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_PASSENGER"], h in HOURS, td in TYPICAL_DAYS}:
	F_t [j, h, td]   = Shares_Mobility_Passenger [j] * (end_uses_input["MOBILITY_PASSENGER"] * mob_pass_time_series [h, td] / t_op [h, td] );
*/	
## Thermal solar & thermal storage:
/*
# [Eq. 26]_V2 relation between dhn thermal solar power and capacity via period capacity factor.
subject to thermal_solar_capacity_factor_dhn {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}, h in HOURS, td in TYPICAL_DAYS}:
	F_t_Solar [j, h, td] <= F_Solar[j] * c_p_t["DHN_SOLAR", h, td];
	
# [Eq. 27]_V2 Overall thermal solar is the sum of specific thermal solar 	
subject to thermal_solar_total_capacity_dhn :
	F ["DHN_SOLAR"] = sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}} (F_Solar[j]);

# [Eq. 28]_V2 Decentralised thermal technology must supply a constant share of heat demand.

subject to decentralised_heating_balance_dhn  {h in HOURS, td in TYPICAL_DAYS} : 
	sum{j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"}} ((F_t [j, h, td] + F_t_Solar [j, h, td]) * layers_in_out[j,"HEAT_LOW_T_DHN"] / cop_time_series[j,"HEAT_LOW_T_DHN", h, td])
	+ sum {i in STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"]} (Storage_out [i,"HEAT_LOW_T_DHN", h, td] - Storage_in [i,"HEAT_LOW_T_DHN", h, td])
	= End_Uses["HEAT_LOW_T_DHN",h,td];*/

/*subject to decentralised_heating_balance_dhn  {h in HOURS, td in TYPICAL_DAYS}:
	sum{j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"] diff {"DHN_SOLAR"},} (F_t [j, h, td] + F_t_Solar [j, h, td]) + sum {i in STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DHN"]} (Storage_out [i, "HEAT_LOW_T_DHN", h, td] - Storage_in [i, "HEAT_LOW_T_DHN", h, td])
		= (pop * end_uses_input["HEAT_LOW_T_HW"] / total_time + (pop * end_uses_input ["HEAT_LOW_T_SH"] ) * heating_time_series [h, td] / t_op [h, td] ) * Share_Heat_Dhn + Network_losses ["HEAT_LOW_T_DHN",h,td] ;
*/

# [Eq. 26] relation between decentralised thermal solar power and capacity via period capacity factor.
subject to thermal_solar_capacity_factor {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"} diff {"DEC_HP_ELEC_A2A"} diff {"DEC_STOVE_WOOD"}, h in HOURS, td in TYPICAL_DAYS}:
	F_t_Solar [j, h, td] <= F_Solar[j] * c_p_t["DEC_SOLAR", h, td];	

# [Eq. 27] Overall thermal solar is the sum of specific thermal solar 	
subject to thermal_solar_total_capacity :
	F ["DEC_SOLAR"] = sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}  diff {"DEC_HP_ELEC_A2A"} diff {"DEC_STOVE_WOOD"}} F_Solar[j];

# [Eq. 28]: Decentralised thermal technology must supply a constant share of heat demand.
subject to decentralised_heating_balance  {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"} diff {"DEC_HP_ELEC_A2A"} diff {"DEC_STOVE_WOOD"} , i in TS_OF_DEC_TECH[j], h in HOURS, td in TYPICAL_DAYS}:
	F_t [j, h, td] + F_t_Solar [j, h, td] + sum {l in LAYERS /*, i in STORAGE_OF_END_USES_TYPES ["HEAT_LOW_T_DECEN"]*/} (Storage_out [i, l, h, td] - Storage_in [i, l, h, td])  
		= Shares_LowT_Dec[j] * (pop * end_uses_input["HEAT_LOW_T_HW"] / total_time + (pop * end_uses_input["HEAT_LOW_T_SH"]) * heating_time_series [h, td] / t_op [h, td]);


## EV storage :
/*
# [Eq. 32] Compute the equivalent size of V2G batteries based on the share of V2G, the amount of cars and the battery capacity per EVs technology
subject to EV_storage_size {j in V2G, i in EVs_BATT_OF_V2G[j]}:
	F [i] = n_car_max * Shares_Mobility_Passenger[j] * Batt_per_Car[j];# Battery size proportional to the amount of cars
	
# [Eq. 33]  Impose EVs to be supplied by their battery.
subject to EV_storage_for_V2G_demand {j in V2G, i in EVs_BATT_OF_V2G[j], h in HOURS, td in TYPICAL_DAYS}:
	Storage_out [i,"ELECTRICITY",h,td] >=  - layers_in_out[j,"ELECTRICITY"]* F_t [j, h, td];*/
		

## Peak demand :

# [Eq. 34] Peak in decentralized heating
subject to peak_lowT_dec {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}, h in HOURS, td in TYPICAL_DAYS}:
	F [j] >= peak_sh_factor * peak_real_factor * worst_year_factor * F_t [j, h, td];

# [Eq. 35] Calculation of max heat demand in DHN (1st constrain required to linearised the max function)
var Max_Heat_Demand >= 0;
subject to max_dhn_heat_demand {h in HOURS, td in TYPICAL_DAYS}:
	Max_Heat_Demand >= End_Uses ["HEAT_LOW_T_DHN", h, td];

# Peak in DHN
subject to peak_lowT_dhn:
	sum {j in TECHNOLOGIES_OF_END_USES_TYPE ["HEAT_LOW_T_DHN"], i in STORAGE_OF_END_USES_TYPES["HEAT_LOW_T_DHN"]} (F [j] + F[i]/storage_discharge_time[i]) >= worst_year_factor * peak_sh_factor * Max_Heat_Demand;


# Limited solar production
subject to roof_lim:
	F["PV"] + (F["DEC_SOLAR"] + F["DHN_SOLAR"])/3 <= 12 * pop;#roof_limit * pop; #https://www.solarthermalworld.org/news/solar-thermal-shows-highest-energy-yield-square-metre


## Adaptation for the case study: Constraints needed for the application to Switzerland (not needed in standard LP formulation)
#-----------------------------------------------------------------------------------------------------------------------

# [Eq. 36]  constraint to reduce the GWP subject to Minimum_gwp_reduction :

/*subject to Minimum_GWP_reduction :
	GWP_op <= pop * gwp_limit;
	#TotalGWP 

# [Eq. 37] Minimum share of RE in primary energy supply
/*subject to Minimum_RE_share :
	sum {j in RE_RESOURCES, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} F_t [j, h, td] * t_op [h, td] 
	>=	re_share_primary *
	sum {j in RESOURCES, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} F_t [j, h, td] * t_op [h, td]	;
*/
	

# [Eq. 38] Definition of min/max output of each technology as % of total output in a given layer. 
subject to f_max_perc {eut in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE[eut]}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [j,h,td] * t_op[h,td]) <= fmax_perc [j] * sum {j2 in TECHNOLOGIES_OF_END_USES_TYPE[eut], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [j2, h, td] * t_op[h,td]);
subject to f_min_perc {eut in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE[eut]}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [j,h,td] * t_op[h,td]) >= fmin_perc [j] * sum {j2 in TECHNOLOGIES_OF_END_USES_TYPE[eut], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [j2, h, td] * t_op[h,td]);


# [Eq. 39] Energy efficiency is a fixed cost
/*subject to extra_efficiency:
	F ["EFFICIENCY"] = 1 / (1 + i_rate);*/


# [Urban] Import when Feed_in
subject to feed_in_tariff:
	sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ( F_t ["ELECTRICITY_FEED_IN", h, td] * t_op [h, td] ) <= sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} ( F_t ["ELEC_EXPORT", h, td] * t_op [h, td] );

# [PoC] Building renovation efficiency upgrade
# subject to building_renovations:
# 	end_uses_input_post["HEAT_LOW_T_SH"] = end_uses_input["HEAT_LOW_T_SH"]*(F ["ROOF"] + F ["FACADE"] + F ["FLOOR"] + F ["WINDOW"]) ; #sum {i in RENOVATION} (F [i])

# [PoC] Definition of min/max output of each technology as % of total output in a given layer. 
subject to renovation_f_max_perc {k in RENOVATION}:
	F [k] <= fmax_perc [k] * end_uses_input ["HEAT_LOW_T_SH"] * pop;

/*
subject to minimum_auto_consumption_rate :
sum{j in HOME_TECHNOLOGIES["ELECTRICITY"], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (layers_in_out[j,"ELECTRICITY"] * F_t [j, h, td] - F_t ["ELEC_EXPORT", h, td])
	>=	auto_consumption_rate*
sum{j in HOME_TECHNOLOGIES["ELECTRICITY"], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (layers_in_out[j,"ELECTRICITY"] * F_t [j, h, td]);
*/
#subject to minimum_auto_sufficiancy_rate :
#sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (End_Uses ["ELECTRICITY", h, td] + End_Uses ["HEAT_LOW_T_DHN", h, td] - F_t ["ELECTRICITY_FEED_IN", h, td] * t_op [h, td] - F_t ["ELECTRICITY", h, td] * t_op [h, td])
#	=	auto_sufficiancy_rate/100*
#sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (End_Uses ["ELECTRICITY", h, td] + End_Uses ["HEAT_LOW_T_DHN", h, td]);
# [PoC] AutoConsumption_Rate if auto_consumption_rate >= 0.3774;  TO BE DEFINED WITH THE INSTALLED CAPACITY
subject to prosumer_policy: 
	Prosumer_tax = (F ["PV"] * 0.085 * 910); # with net-metering ( à ajouter dans totalcosts) #en considérant que la capa installée est en kwe ou = kwc ??



##########################
### OBJECTIVE FUNCTION ###
##########################

# Can choose between TotalGWP and TotalCost
minimize obj: TotalCost;

/*

##############################################
###            GLPK version                ###
##############################################

solve;

*/

### Printing output



