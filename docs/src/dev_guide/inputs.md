# Inputs

Inputs used to run the holistic electricity model: [indices](#indices) and [parameters](#parameters). [Indices](#indices) refer to the dimensions using which the rest of the data is organized, such as the number of generation and storage resources, independent power producers, the balancing areas condsidered, the representative days considered etc.

## Indices

### Common Indices

 - **index_y.csv**: Representative simulation years. It consisits of a single row of year values. The file is used to generate two types of indices, `index_y` which is the years for a simulation window and `index_y_fix` which represents the full simulation horizon (does not change).

 - **index_s.csv**: New resource depreciation years represented as a row vector.

 - **index_d.csv**: Representative day and hour (from ReEDS). Each day is represented as `d[0-9][0-9][0-9]` in a row vector.

 - **index_t.csv**: Representative hour within each representative day (from ReEDS) represented as `h[0-9][0-9][0-9]` in a row vector.

 - **index_h.csv**: Customer groups. Example values are `Residential, Commercial, Industrial`.

 - **index_j.csv**: Green technology types. Example values are `biopower, wind-ons, wind-ofs`.

 - **index_z.csv**: Balancing Area zones as modeled in ReEDS. Refer to the [model regions map](https://www.nrel.gov/analysis/assets/pdfs/reeds-model-regions-map.pdf).

 - **index_sector.csv**: Customer high level groups. Example values are `Residential, Commercial, Industrial`.

 - **index_h_sector_mapping.csv**: File that maps each zone in index\_z to a group in index\_sector. It conisits of two columns, `index_z` and `index_h`.

### Agent-Specfic Indices

#### Regulator

 - **index_rate_tou.csv**: Index for supported time-of-use rates. Example values are `Off-Peak, On-Peak, Mid-Peak`.

 - **tou_rate_structure_\$(tou_suffix).csv**: Time-of-use rate for each representative hour and representative day combination. The column names are `index_d`, `index_t` and `index_rate_tou`. Here, the input file selected depends on the `tou_suffix` specified in the config.yaml file under `regulator_options`.

#### VIU/IPP

 - **index_k_existing.csv**: Existing bulk generation technologies. Consists of a single row containing the existing technolgies, for e.g., `coaloldscr`, `gas-ct`, `nuclear`, `wind-ofs` etc.

 - **index_k_new.csv**: Potential bulk generation technologies. Consists of a single row containing the new technolgies, for e.g., `coaloldscr`, `gas-ct`, `nuclear`, `wind-ofs` etc.

 - **index_rps.csv**: [Renewable Potfolio Standard (RPS)](https://www.nrel.gov/state-local-tribal/basics-portfolio-standards.html) qualified technologies. Consists of a single row containing qualified technolgies, for e.g., `biopower`, `wind-ofs` etc.

 - **index_stor_existing.csv**: Existing bulk storage technologies. Consists of a single row containing the existing technolgies, for e.g., `battery_1`, `battery_2` etc.

 - **index_stor_new.csv**: Potential bulk storage technologies. Consists of a single row containing the existing technolgies, for e.g., `battery_1`, `battery_2` etc.

 - **index_l.csv**: Transmission lines. Each line is represented as `\$(origin_zone)||\$(destination_zone)`, for e.g., `p100||p101`.

 - **index_p.csv**: Individual independent power producers (IPP). Cosnsists of a single row of IPP names.

#### Customer

 - **index_m.csv**: Behind-The-Meter (BTM) technologies. Eaxample values are `BTMPV, BTMStorage`.

## Parameters

### Model Parameters

 - **Omega.csv**: Number of days per representative day. Column names are the representative days in index_d and the values are the number of days for that representative day.

 - **Gamma.csv**: Number of customers of type in index\_h at zone in index\_z. Column names are index\_z, followed by the customer groups (for e.g., `index_z, Commercial, Residential, Industrial`).

 - **Year.csv**: Column names and values are the representative years in index_y.

 - **Time.csv**: Number of days per representative hour. Column names are the representative hours in index_t and the values are the number of hours for that representative hour.

 ### Hourly Resource Data

 - **AvailabilityDER.csv**: Availability of DER units (fraction) by customer type, BTM technology, zone, representative day and hour. Column names are index\_h, index\_m, index\_z, index\_d, followed by represenative hours in index\_t  (for e.g., `index_h, index_m, index_z, index_d, h001, h002` etc.).

 - **AvailabilityNewIPP.csv**: Availability of new units (%) by IPP type, technology type, zone, representative day and hour. Column names are index\_p, index\_k\_new, index\_z, index\_d, followed by represenative hours in index\_t  (for e.g., `index_p, index_k_existing, index_z, index_d, h001, h002` etc.).

 - **AvailabilityOldIPP.csc**: Availability of existing units (fraction) by IPP type, technology type, zone, representative day and hour. Column names are index\_p, index\_k\_existing, index\_z, index\_d, followed by represenative hours in index\_t  (for e.g., `index_p, index_k_existing, index_z, index_d, h001, h002` etc.).


### Financial Parameters

#### New Resources

 - **Lifetime_new.csv**: Lifetime of new generation technologies in years. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **LifetimeNewIPP.csv**: Lifetime of new generation technologies in years by IPPs in index\_p. Columns are index\_p, followed by technology names (for e.g., `index_p, gas-ct, gas-cc, upv` etc.).

 - **LifetimeStorNewIPP.csv**: Lifetime of new storage technologies in years by IPPs in index\_p. Columns are index\_p, followed by storage technology names (for e.g., `index_p, battery_1, battery_2` etc.).

#### Incentives

 - **PTCNew.csv**: Production Tax Credit (PTC) for new generation technologies in \$/MW. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **PTCNewmy.csv**: Multi-year Production Tax Credit (PTC) for new generation technologies in \$/MW. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **PTCOld.csv**: Production Tax Credit (PTC) for existing generation technologies in \$/MW. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **PTCNewmy.csv**: Multi-year Production Tax Credit (PTC) for new generation technologies in \$/MW. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **ITCNew.csv**: Investment Tax Credit (ITC) (fraction) for new generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **ITCNewmy.csv**: Multi-year Investment Tax Credit (ITC) (fraction) for new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **ITCStorNewmy.csv**: Multi-year Investment Tax Credit (ITC) (fraction) for new storage technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, battery_1, battery_2` etc.).

 - **ITCOld.csv**: Investment Tax Credit (ITC) (fraction) for existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **ITCStorOld.csv**: Investment Tax Credit (ITC) (fraction) for existing storage technologies. Columns are technology names (for e.g., `battery_1, battery_2` etc.).

#### Bulk(VIU/ IPP) Resources

 - **DebtRatio.csv**: Debt ratio. Currently it is only set from a file for IPP (for others they are constants). Column name is the IPP agent name (for e.g. `ipp1`).

 - **COD.csv**: Cost of debt (%). Currently it is only used for IPP (for others they are constants). Column name is the IPP agent names in index\_p (for e.g. `ipp1`).

 - **COE.csv**: Cost of equity (fraction). Currently it is only used for IPP (for others they are constants). Column name is the IPP agent names in index\_p (for e.g. `ipp1`).

 - **Tax.csv**: Tax rate by IPP type. Column names are the IPP names in index\_p (for e.g., `ipp1`).

 - **CumuTaxDepreOld.csv**: Cumulative tax depreciation of existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuTaxDepreOldmy.csv**: Multi-year cumulative tax depreciation of existing generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **CumuTaxDepreStorOldmy.csv**: Multi-year cumulative tax depreciation of existing storage technologies. Columns are index\_y, followed by existing storage technology names (for e.g., `index_y, battery_1, battery_2` etc.).

 - **CumuTaxDepreNew.csv**: Cumulative tax depreciation of new generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuTaxDepreNewmy.csv**: Multi-year cumulative tax depreciation of new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **CumuTaxDepreStorNewmy.csv**: Multi-year cumulative tax depreciation of new storage technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, battery_1, battery_2` etc.).

 - **CumuAccoutDepreOld.csv**: Cumulative tax depreciation of existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuAccoutDepreNew.csv**: Cumulative accounting depreciation of new generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuAccoutDepreNewmy.csv**: Multi-year accounting depreciation of new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **CumuAccoutDepreStorNewmy.csv**: Multi-year accounting depreciation of new storage technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **AnnualTaxDepreNewmy.csv**: Multi-year annual tax depreciation of new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **AnnualTaxDepreStorNewmy.csv**: Multi-year annual tax depreciation of new storage technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).
 
 - **AnnualTaxDepreOldmy.csv**: Multi-year annual tax depreciation of existing generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **AnnualTaxDepreStorOldmy.csv**: Multi-year annual tax depreciation of existing storage technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, battery_1, battery_2` etc.).

 - **AnnualAccoutDepreNewmy.csv**: Multi-year annual accounting depreciation of existing generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **AnnualAccoutDepreStorNewmy.csv**: Multi-year annual accounting depreciation of new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, battery_1, battery_2` etc.).

 - **AnnualAccoutDepreOld.csv**: Annual accounting depreciation of existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuAccoutDepreOldmy.csv**: Multi-year accounting depreciation of existing generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **CumuAccoutDepreStorOldmy.csv**: Multi-year accounting depreciation of existing storage technologies. Columns are index\_y, followed by storage technology names (for e.g., `index_y, battery_1, battery_2` etc.).

 - **CumuITCAmortNew.csv**: Cumulative amortized ITC (fraction) for new generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuITCAmortNewmy.csv**: Multi-year cumulative amortized ITC (fraction) for new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **CumuITCAmortStorNewmy.csv**: Multi-year cumulative amortized ITC (fraction) for new generation technologies. Columns are index\_y, followed by technology names (for e.g., `index_y, gas-ct, gas-cc, upv` etc.).

 - **CumuITCAmortOldmy.csv**: Cumulative amortized ITC (fraction) for existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **CumuITCAmortStorOldmy.csv**: Cumulative amortized ITC (fraction) for existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **AnnualITCAmortNewmy.csv**: Annual amortized ITC (fraction) for new generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **AnnualITCAmortStorNewmy.csv**: Annual amortized ITC (fraction) for new generation technologies. Columns are technology names (for e.g., `battery_1, battery_2` etc.).

 - **AnnualITCAmortOldmy.csv**: Annual amortized ITC (fraction) for existing generation technologies. Columns are technology names (for e.g., `gas-ct, gas-cc, upv` etc.).

 - **AnnualITCAmortStorOldmy.csv**: Annual amortized ITC (fraction) for existing generation technologies. Columns are technology names (for e.g., `battery_1, battery_2` etc.).


#### Demand-Side (Customer Resources)


 - **CustomerDebtRatio.csv**: Customer debt ratio by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (for e.g., `index_z, Residential, Commercial, Industrial`).

 - **CustomerCOD.csv**: Customer cost of debt (fraction) by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (for e.g., `index_z, Residential, Commercial, Industrial`).

 - **CustomerCOE.csv**: Customer cost of equity (fraction) by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (for e.g., `index_z, Residential, Commercial, Industrial`).

 - **CustomerTax.csv**: Customer tax rate by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (for e.g., `index_z, Residential, Commercial, Industrial`).

 - **DER_ITCNew.csv**: ITC for new DERs. Columns are BTM technology names in index\_m (for e.g., `BTMPV, BTMStorage`).

 - **DER_ITCNewmy.csv**: Multi-year ITC for new DERs. Columns are index\-y, followed by BTM technology names in index\_m (for e.g., `index_y, BTMPV, BTMStorage`).

#### SupplyChoice Use Case and GreenDeveloper (not used)

 - **RetailCompetition.csv**: Retail competition rates by year (used to calculate green tariff subscription percentages). Column names are years in index\_y (for e.g., `2021, 2022` etc.).

#### Rate-Making/ Regulator Costs

 - **distribution_cost.csv**: Distribution cost in \$ for each index\_z zone and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (For e.g., `index_y, 2021, 2022`). The filenames imply the corresponding type of cost.

 - **administration_cost.csv**: Administration cost in \$ for each index\_z zone and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (For e.g., `index_y, 2021, 2022`). The filenames imply the corresponding type of cost.

 - **transmission_cost.csv**: Transmission cost in \$ for each index\_z zone and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (For e.g., `index_y, 2021, 2022`). The filenames imply the corresponding type of cost.

 - **interconnection_cost.csv**: Interconnection cost in \$ for each index\_z zone and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (For e.g., `index_y, 2021, 2022`). The filenames imply the corresponding type of cost.

 - **system_cost.csv**: System cost in \$ for each index\_z zone and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (For e.g., `index_y, 2021, 2022`). The filenames imply the corresponding type of cost.

#### Bulk System (VIU/IPP) Cost Parameters

 - **FOM_new.csv**: Fixed Operation and Maintenance cost of new generation technologies in \$/MW-yr for each index\_z zone. Columns are index\_z, followed by technology names (for e.g., `index_z, gas-ct, gas-cc, upv` etc.).

 - **FOMNewIPP.csv**: Fixed Operation and Maintenance cost of new generation technologies in \$/MW-yr for each IPP in index\_p. Columns are index\_p, followed by technology names (for e.g., `index_p, gas-ct, gas-cc, upv` etc.).

 - **FOMNewIPPmy.csv**: Multi-year fixed Operation and Maintenance cost of new generation technologies in \$/MW-yr for each IPP in index\_p. Columns are index\_y, index\_p, followed by technology names (for e.g., `index_y, index_p, gas-ct, gas-cc, upv` etc.).

 - **CapEx_new.csv**: Capital cost of new generation technologies in \$/MW for each index\_z zone. Columns are index\_z, followed by technology names (for e.g., `index_z, gas-ct, gas-cc, upv` etc.).

 - **CapExNewIPP.csv**: Capital cost of new generation technologies in \$/MW for each ipp in index\_p. Columns are index\_p, followed by technology names (for e.g., `index_p, gas-ct, gas-cc, upv` etc.).

 - **CapExNewIPPmy.csv**: Multi-year capital cost of new generation technologies in \$/MW for each ipp in index\_p and year in index\_y. Columns are index\_y, index\_p, followed by technology names (for e.g., `index_y, index_p, gas-ct, gas-cc, upv` etc.).

 - **StorCapExNewIPPmy.csv**: Multi-year capital cost of new storage technologies in \$/MW for each ipp in index\_p and year in index\_y. Columns are index\_y, index\_p, followed by technology names (for e.g., `index_y, index_p, battery_1, battery_2` etc.).

 - **CapExOld.csv**: Capital cost of existing generation technologies in \$/MW for each index\_z zone. Columns are index\_z, followed by technology names (for e.g., `index_z, gas-ct, gas-cc, upv`etc.).

 - **CapExStorOld.csv**: Capital cost of existing storage technologies in \$/MW for each index\_z zone. Columns are index\_z, followed by storage technologies (for e.g., `index_z, battery_1, battery_2`).

 - **FixedCostOld.csv**: Fixed cost of existing generation technologies in \$/MW by zone. Column names are index\_z, followed by technology names (for e.g., `index_z, gas-ct, gas-cc, upv` etc.).

 - **FixedCostOldIPPmy.csv**: Multi-year fixed cost of existing generation technologies in \$/MW by zone. Column names are index\_y, index\_z, followed by technology names (for e.g., `index_y, index_p, gas-ct, gas-cc, upv` etc.).

 - **FixedCostStorOldmy.csv**: Multi-year fixed cost of existing storage technologies in \$/MW by zone. Column names are index\_y, index\_z, followed by technology names (for e.g., `index_y, index_z, battery_1, battery_2` etc.).

 - **FixedCostOldIPP.csv**: Fixed cost of existing generation technologies in \$/MW by IPP and zone. Column names are index\_p, index\_z, followed by technology names (for e.g., `index_p, index_z,gas-ct, gas-cc, upv` etc.).

 - **FixedCostOldIPPmy.csv**: Multi-year fixed cost of existing generation technologies in \$/MW by IPP and zone. Column names are index\_y, index\_p, index\_z, followed by technology names (for e.g., `index_y, index_p, index_z,gas-ct, gas-cc, upv` etc.).

 - **FixedCostStorOldIPPmy.csv**: Multi-year fixed cost of existing storage technologies in \$/MW by IPP and zone. Column names are index\_y, index\_p, index\_z, followed by technology names (for e.g., `index_y, index_p, index_z,battery_1, battery_2` etc.).

 - **VariableCostNew.csv**: Variable cost of new generation technologies in \$/MW by zone. Column names are index\_z, followed by technology names (for e.g., `index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostNewmy.csv**: Multi-year Variable cost of new generation technologies in \$/MW by zone. Column names are index\_y, index\_z, followed by technology names (for e.g., `index_y, index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostOld.csv**: Variable cost of existing generation technologies in \$/MW by zone. Column names are index\_z, followed by technology names (for e.g., `index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostOldmy.csv**: Multi-year Variable cost of existing generation technologies in \$/MW by zone. Column names are index\_y, index\_z, followed by technology names (for e.g., `index_y, index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostNewIPP.csv**: Variable cost of new generation technologies in \$/MW by IPP and zone. Column names are index\_p, index\_z, followed by technology names (for e.g., `index_p, index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostNewIPPmy.csv**: Multi-year Variable cost of new generation technologies in \$/MW by IPP and zone. Column names are index\_y, index\_p, index\_z, followed by technology names (for e.g., `index_y, index_p, index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostOldIPP.csv**: Variable cost of existing generation technologies in \$/MW by IPP and zone. Column names are index\_p, index\_z, followed by technology names (for e.g., `index_p, index_z,gas-ct, gas-cc, upv` etc.).

 - **VariableCostOldIPPmy.csv**: Multi-year Variable cost of existing generation technologies in \$/MW by IPP and zone. Column names are index\_y, index\_p, index\_z, followed by technology names (for e.g., `index_y, index_p, index_z,gas-ct, gas-cc, upv` etc.).

 - **StorFOMNewIPPmy.csv**: Multi-year fixed Operation and Maintenance cost of new storage technologies in \$/MW-yr for each IPP in index\_p. Columns are index\_y, index\_p, followed by technology names (for e.g., `index_y, index_p, battery_1, battery_2` etc.).

#### Customer Cost Parameters

 - **CapExDER.csv**: Capital cost of BTM technologies in \$/MW for each index\_z zone and customer type in index\_h. Columns are index\_z, index\_h, followed by technologies (for e.g., `index_z, index_h, BTMPV, BTMStorage`).

 - **CapExDERmy.csv**: Multi-year capital cost of BTM technologies in \$/MW for each year in index\_y, index\_z zone and customer type in index\_h. Columns are index\_y, index\_z, index\_h, followed by technologies (for e.g., `index_y, index_z, index_h, BTMPV, BTMStorage`).

 - **FOMDER.csv**: Fixed Operation and Maintenance cost of new BTM generation technologies in \$/MW-yr for each zone in index\_z and for each customer type in index\_h. Columns are index\_z, index\_h, followed by technology names (for e.g., `index_z, index_h, BTMPV, BTMStorage` etc.).

#### Cost Parameter applicable to WM/ IPP + ISO for the Capacity Market

 - **NetCONE.csv**: Net Cost of New Entry for capacity market in in \$/MW-yr. Column names are the years (for e.g., `2021, 2022` etc.).

### Technology Parameters

#### Technology Performance Assumptions for Existing and Future Resources

 - **StorRTENewIPPmy.csv**: Multi-year round-trip efficiency (fraction) of new storage resources by year, IPP and zone. Column names are index\_y, index\_p, index\_z, followed by the technolgies (for e.g., `index_y, index_p, index_z, battery_1, battery_2` etc.).
 
 - **StorRTEOldIPPmy.csv**: Multi-year round-trip efficiency (fraction) of exisiting storage resources by year, IPP and zone. Column names are index\_y, index\_p, index\_z, followed by the technolgies (for e.g., `index_y, index_p, index_z, battery_1, battery_2` etc.).

 - **NewStorDuration.csv**: Duration in hours of new storage technologies. Column names are the new storage technology names (for e.g., `battery_1, battery_2` etc.).

#### Technology Performance Assumptions for Customers

 - **rte_dist_stor.csv**: Distributed storage round trip efficiency (fraction) for a customer type in a zone. Column names are index\_z, followed by the customer type values in index\_h (for e.g., `index_z, Commercial, Residential, Industrial`).

 - **duration_dist_stor.csv**: Distributed storage duration for a customer type in a zone in hours. Column names are index\_z, followed by the customer type values in index\_h (for e.g., `index_z, Commercial, Residential, Industrial`).

#### Technology Performance for Existing and Future Resources

 - **CapacityCredit_new.csv**: Capacity credit (fraction) of existing generation technologies by year and zone. Column names are index\_y, index\_z, followed by the existing technology names (for e.g., `index_y, index_z, gas-ct, gas-cc, upv` etc.).

 - **CapacityCreditStor_new.csv**: Capacity credit (fraction) of new storage technologies by year and zone. Column names are index\_y, index\_z, followed by the new storage technology names (for e.g., `index_y, index_z, battery_1, battery_2` etc.).

 - **CapacityCredit_old.csv**: Capacity credit (fraction) of existing generation technologies by year and zone. Column names are index\_y, index\_z, followed by the existing technology names (for e.g., `index_y, index_z, gas-ct, gas-cc, upv` etc.).

 - **CapacityCreditStor_old.csv**: Capacity credit (fraction) of existing storage technologies by year and zone. Column names are index\_y, index\_z, followed by the existing storage technology names (for e.g., `index_y, index_z, battery_1, battery_2` etc.).


 - **EmissionRateNewIPPmy.csv**: Multi-year emission rate in pounds per million Btu by year, IPP, zone and new technology types. Column names are index\_y, index\_p, index\_z, followed by technolgies (for e.g., `index_y, index_p, index_z, gas-ct, gas-cc, upv` etc.).

 - **EmissionRateOldIPPmy.csv**: Multi-year emission rate in pounds per million Btu by year, IPP, zone and existing technology types. Column names are index\_y, index\_p, index\_z, followed by technolgies (for e.g., `index_y, index_p, index_z, gas-ct, gas-cc, upv` etc.).

- **ExistingStorInitialEnergyIPP.csv**: Initial energy in MW for each IPP, zone and existing storage generation type. Column names are index\_p, index\_z, followed by existing storage technologies (for e.g., `index_p, index_z, battery_1, battery_2` etc.).

- **NewStorInitialEnergyIPP.csv**: Initial energy in MW for each IPP, zone and new storage generation type. Column names are index\_p, index\_z, followed by new storage technologies (for e.g., `index_p, index_z, battery_1, battery_2` etc.).

#### Technology Performance for Customers

 - **OptimalDER.csv**: Optimal DER capacity in MWh of the BTM technologies for each customer group in each zone. Column names are index\_z, index\_h, folowed by the technology names in index\_m (for e.g., `index_z, index_h, BTMPV, BTMStorage`).

 - **OptimalDERmy.csv**: Multi-year optimal DER capacity in MWh of the BTM technologies for each customer group in each zone in each year. Column names are index\_y, index\_z, index\_h, folowed by the technology names in index\_m (for e.g., `index_y, index_z, index_h, BTMPV, BTMStorage`).

 - **initial\_energy\_dist\_stor.csv**: Initial ditributed storage energy in MWh for a customer type in a zone. Column names are index\_z, followed by the customer type values in index\_h (for e.g., `index_z, Commercial, Residential, Industrial`).

#### Technology Adoption Parameters (Specific to Customers)

 - **Bass_P.csv**: Coefficient of innovation (p) values for the Bass diffusion model (used to construct the PV adoption model) for each cutomer type, zone and BTM technology combination. Column names are index\_z, index\_h, followed by the values in index\_m (for e.g., `index_z, index_h, BTMPV, BTMStorage`).

 - **Bass_Q.csv**: Coefficient of imitation (q) values for the Bass diffusion model (used to construct the PV adoption model) for each cutomer type, zone and BTM technology combination. Column names are index\_z, index\_h, followed by the values in index\_m (for e.g., `index_z, index_h, BTMPV, BTMStorage`).

 - **DC_length.csv**: Capacity market demand curve length for each year in index\_y. Column names are years (for e.g., `2021, 2022`).

 - **WTP.csv**: Wllingness to pay in \$/MW for each year in index\_y. Column names are the years in index\_y (for e.g., `2021, 2022` etc.).

### Existing System Parameters

 - **ExistingStorDuration.csv**: Duration in hours of existing storage technologies. Column names are the existing storage technology names (for e.g., `battery_1, battery_2` etc.).

 - **TransmissionTopology.csv**: Topology of the transmission system expressed as an incidence matrix between the lines and the zones. Column names are index\_l, followed by the zones in index\_z (for e.g., `index_l, p100, p101` etc.) with the values being +1 if the zone is the origin or -1 if it is the destination.

 - **ExistingCapacityIPP.csv**: Existing capacity in MW for each IPP, zone and generation type. Column names are index\_p, index\_z, followed by existing generation technologies (for e.g., `index_p, index_z, gas-ct, gas-cc, upv` etc.).

 - **ExistingStorCapacityIPP.csv**: Storage capacity in MW for each IPP, zone and existing generation type. Column names are index\_p, index\_z, followed by existing storage technologies (for e.g., `index_p, index_z, battery_1, battery_2` etc.).

 - **ExistingDER.csv**: Exisiting DER capacity in MWh of the BTM technologies for each customer group in each zone. Column names are index\_h, index\_z, folowed by the technology names in index\_m (for e.g., `index_h, index_z, BTMPV, BTMStorage`).

 - **Export.csv**: Net export in MWh in the representative hours for each index\_z zone and index\_d day combination. The column names are index\_z, index\_d, followed by the representative hours (for e.g., `index_z, index_d, h001, h002, h003`).

 - **Exportmy.csv**: Multi-year net export in MWh in the representative hours for each index\_y year, index\_z zone, index\_d day combination. The column names are index\_y, index\_z, index\_d, followed by the representative hours (for e.g., `index_y, index_z, index_d, h001, h002, h003`).

 - **Demand.csv**: Demand in MWh for each customer group, zone and representative day in each representative hour. The column names are index\_h, index\_z, index\_d, followed by the representative hours (for e.g., `index_h, index_z, index_d, h001, h002, h003`).

 - **Demandmy.csv**: Multi-year demand in MWh for each customer group, zone, year and representative day in each representative hour. The column names are index\_h, index\_h, index\_z, index\_d, followed by the representative hours (for e.g., `index_y, index_h, index_z, index_d, h001, h002, h003`).