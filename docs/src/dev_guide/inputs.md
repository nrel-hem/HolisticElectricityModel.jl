# Inputs

Inputs used to run the Holistic Electricity Model: [indices](#indices) and [parameters](#parameters).

## Indices
Indices are the dimensions along which the rest of the data are organized, such as the types of generation and storage resources, number of independent power producers, the balancing areas considered, the representative days and hours in the model, etc. Indices are always single-dimensional, with individual values arrayed in the first row of a .csv file (i.e., serialized as a row vector).

In limited cases, mappings between two or three indices must also be provided. Mapping files contain two or three columns labeled with the name of each index being mapped. The presence of an entry is interpreted as membership or a scaling factor of 1.0, whereas absence from the file is interpreted as non-membership or a scaling factor of 0.0.

### Common Indices

 - **index_y.csv**: Simulation years, e.g., 2020, 2021, .... The file is used to generate two indices within HEM. The current simulation window is defined by `index_y`. `index_y_fix` represents the full simulation horizon (equal to the years listed in the .csv file) and does not change. Year indices are used in HEM as both strings (for labeling other data) and integers (for iterating and performing calculations).

 - **index_s.csv**: Years from initial build over which new resources depreciate. For example, this could be a list of integers 1, 2, ..., 15.

 - **index_d.csv**: Representative days. When using data from ReEDS, each day is represented as `d[0-9][0-9][0-9]` but there is no constraint in the model regarding how representative days are labeled.

 - **index_t.csv**: Row vector of representative hours within each representative day. When using data from ReEDS, each hour is represented as `h[0-9][0-9][0-9]`. Currently the model assumes that all characters except the first in each label can be converted to integer and interpreted as the last hour of the day (of the hours 1 to 24) included in this representative time period. The model also assumes that each representative hour represents the same amount of time.

 - **index_h.csv**: Customer groups. Example values are `Residential, Commercial, Industrial`, but this index can also be used to differentiate customer data more finely.

 - **index_z.csv**: Balancing areas (BAs) / utility zones. [HEMData.jl](https://github.com/nrel-hem/HolisticElectricityModelData.jl) provides data parsing tied to the [ReEDS BAs]((https://www.nrel.gov/analysis/assets/pdfs/reeds-model-regions-map.pdf).

 - **index_sector.csv**: High-level customer groups used to establish classes for rate-making. Example values are `Residential, Commercial, Industrial`.

 - **index_h_sector_mapping.csv**: File that maps each customer type in index\_h to a high-level customer group (for rate-making) in index\_sector. It consists of two columns, `index_h` and `index_sector`.

### Agent-Specfic Indices

#### Regulator

 - **index_rate_tou.csv**: Index for supported time-of-use rate periods. Example values are `Off-Peak, On-Peak, Mid-Peak`.

 - **tou_rate_structure_\$(tou_suffix).csv**: Time-of-use rate period assigned to each representative day and hour combination. The column names are `index_d`, `index_t` and `index_rate_tou`. Here, the input file selected depends on the `tou_suffix` specified in the config.yaml file under `regulator_options`.

#### Utility and IPPGroup

 - **index_k_existing.csv**: Existing bulk generation technologies, e.g., `coaloldscr`, `gas-ct`, `nuclear`, `wind-ofs` etc.

 - **index_k_new.csv**: Bulk generation technologies available for new investment, e.g., `coaloldscr`, `gas-ct`, `nuclear`, `wind-ofs` etc.

 - **index_rps.csv**: [Renewable Potfolio Standard (RPS)](https://www.nrel.gov/state-local-tribal/basics-portfolio-standards.html) qualified technologies, e.g., `biopower`, `wind-ofs` etc.

 - **index_stor_existing.csv**: Existing bulk storage technologies, e.g., `battery_1`, `battery_2` etc.

 - **index_stor_new.csv**: Bulk storage technologies available for new investment, e.g., `battery_1`, `battery_2` etc.

 - **index_l.csv**: Transmission lines. Each line is represented as `\$(origin_zone)||\$(destination_zone)`, e.g., `p100||p101`. Both `origin_zone` and `destination_zone` are required to be elements of `index_z`.

 - **index_p.csv**: Individual independent power producers (IPP). Cosnsists of a single row of IPP names.

#### Customer

 - **index_m.csv**: Behind-The-Meter (BTM) technologies, e.g., `BTMPV, BTMStorage`. Currently these index values are hard-coded into the model and would be difficult to change/remove.

#### Green Developer

 - **index_j.csv**: Green technology types. Example values are `biopower, wind-ons, wind-ofs`. This index is related to the `GreenDeveloper` agent and direct purchases of clean electricity by customers through power purchase agreements.

## Parameters

Parameters contain numeric data defined over zero, one, or more indices. If the parameters is a scalar (no indices), its value is expected to be in the first row and column of the .csv file. If the parameter is defined over a single index, the elements of that index (values or tags that define the contents of the index) are listed in the first row, with the corresponding numeric values listed below, in the second row of the .csv file. If the parameter is defined over multiple (`n`) indices, all but one of the indices form the first `n-1` columns, with the name of the index in the first row and its elements listed below. The last index's elements are listed in the first row starting in the $n^{th}$ column. Numeric values are arrayed under the columns defined by the last index's elements.

### Model Parameters

 - **Omega.csv**: Number of days per representative day. Defined over index\_d, and generally expected to sum to 365.0.

 - **Gamma.csv**: Number of customers of type h (index\_h) in zone z (index\_z). Column names are index\_z, followed by the index\_h elements (e.g., `Commercial, Residential, Industrial`).

 - **Year.csv**: Numeric values of index\_y, defined over index\_y, that is, two rows containing the elements in index\_y.

 - **Time.csv**: Last hour (of hours 1, 2, ..., 24) represented by each element of index\_t, defined over index\_t. For ReEDS representative hours, corresponds to the integer obtained by parsing all but the first character of the element string (i.e., 4, 8, 12, 16, 20, 24).


### Existing System Parameters (used throughout, especially by Utility, IPPGroup, CustomerGroup, and Regulator)

 - **ExistingStorDuration.csv**: Duration (hours) of new storage technologies. Column names are the new storage technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **ExistingStorInitialEnergy.csv**: Initial stored energy (Mwh) for existing storage resources in index\_stor\_existing by year in index\_y, zone in index\_z and day in index\_d. Column names are index\_y, index\_stor, index\_z, followed by the days in index\_d (e.g., `d001, d002, ...`).

 - **ExistingStorInitialEnergyIPP.csv**: Initial energy (MWh) for existing storage resources in index\_stor\_existing by year in index\_y, IPP in index\_p, zone in index\_z and day in index\_d. Column names are index\_y, index\_stor, index\_p, index\_z, followed by the days in index\_d (e.g., `d001, d002, ...`).

 - **TransmissionTopology.csv**: Topology of the transmission system expressed as an incidence matrix between the lines in index\_l and the zones in index\_z. Column names are index\_l, followed by the zones in index\_z (e.g., `p100, p101, ...` etc.) with the values being +1 if the zone is the origin or -1 if it is the destination.

 - **ExistingCapacityIPP.csv**: Existing capacity (MWh) for each IPP in index\_p, zone in index\_z and generation type in index\_k\_existing. Column names are index\_p, index\_z, followed by existing generation technologies in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **ExistingStorCapacityIPP.csv**: Storage capacity (MWh) for each IPP in index\_p, zone in index\_z and generation type in index\_stor\_existing. Column names are index\_p, index\_z, followed by existing storage technologies in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **ExistingDER.csv**: Exisiting DER capacity (MWh) of the BTM technologies for each customer group in index\_h and zone in index\_z. Column names are index\_h, index\_z, folowed by the technology names in index\_m (e.g., `BTMPV, BTMStorage`).

 - **ExistingDERmy.csv**: Exisiting DER capacity (MWh) of the BTM technologies for each customer group in index\_h and zone in index\_z. Column names are index\_h, index\_z, folowed by the technology names in index\_m (e.g., `BTMPV, BTMStorage`).

 - **Export.csv**: Net export (MWh) in the representative hours for each index\_z zone and index\_d day combination. The column names are index\_z, index\_d, followed by the representative hours in index\_t (e.g., `h001, h002, h003`).

 - **Exportmy.csv**: Multi-year net export (MWh) in the representative hours for each index\_y year, index\_z zone, index\_d day combination. The column names are index\_y, index\_z, index\_d, followed by the representative hours in index\_t (e.g., `h001, h002, h003`).

 - **Demand.csv**: Demand (MWh) for each customer group in index\_h, zone in index\_z, representative day in index\_d and representative hour in index\_t. The column names are index\_h, index\_z, index\_d, followed by the representative hours in index\_t (e.g., `h001, h002, h003`).

 - **Demandmy.csv**: Multi-year demand (MWh) for each year in index\_y, customer group in index\_h, zone in index\_z, representative day in index\_d and representative hour in index\_t. The column names are index\_y, index\_h, index\_z, index\_d, followed by the representative hours in index\_t (e.g., `h001, h002, h003`).


### Hourly Resource Data

 - **AvailabilityDER.csv**: Availability of DER units (fraction) defined over customer types in index\_h, BTM technologies in index\_m, zones in index\_z, representative days in index\_d and hours in index\_t. Column names are index\_h, index\_m, index\_z, index\_d, followed by represenative hours in index\_t (e.g., `h004, h008, ...`).

 - **AvailabilityNew.csv**: Availability of new units (fraction) defined over technology types in index\_k\_new, zones in index\_z, representative days in index\_d and hours in index\_t. Column names are index\_k\_new, index\_z, index\_d, followed by represenative hours in index\_t  (e.g., `h001, h002, ...`).

 - **AvailabilityNewIPP.csv**: Availability of new units (fraction) defined over IPPs in index\_p, technology types in index\_k\_new, zones in index\_z, representative days in index\_d and hours in index\_t. Column names are index\_p, index\_k\_new, index\_z, index\_d, followed by represenative hours in index\_t  (e.g., `h001, h002, ...`).

 - **AvailabilityOld.csv**: Availability of existing units (fraction) defined over technology types in index\_k\_existing, zones in index\_z, representative days in index\_d and hours in index\_t. Column names are index\_k\_new, index\_z, index\_d, followed by represenative hours in index\_t  (e.g., `h001, h002, ...`).

 - **AvailabilityOldIPP.csc**: Availability of existing units (fraction) defined over IPPs in index\_p, technology types in index\_k\_existing, zones in index\_z, representative days in index\_d and hours in index\_t. Column names are index\_p, index\_k\_existing, index\_z, index\_d, followed by represenative hours in index\_t  (e.g., `h001, h002, ...`).


### Financial Parameters

#### New Resources

 - **Lifetime_new.csv**: Lifetime of new generation technologies in index\_k\_new in years, columns are technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv` etc.).

 - **LifetimeNewIPP.csv**: Lifetime of new generation technologies in years by IPPs in index\_p. Columns are index\_p, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **LifetimeStorNew.csv**: Lifetime of storage technologies in index\_k\_new in years, columns are technology names in index\_stor\_new (e.g., `battery_2, battery_4, ...`).

 - **LifetimeStorNewIPP.csv**: Lifetime of new storage technologies in years by IPPs in index\_p. Columns are index\_p, followed by storage technology names in index\_stor\_new(e.g., `battery_1, battery_2, ...`).

#### Incentives

 - **PTCNew.csv**: Production Tax Credit (PTC) (\$/MW) for new generation technologies. Columns are technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **PTCNewmy.csv**: Multi-year Production Tax Credit (PTC) (\$/MW) for new generation technologies defined over years in in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **PTCOld.csv**: Production Tax Credit (PTC) (\$/MW) for existing generation technologies. Columns are technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **PTCNewmy.csv**: Multi-year Production Tax Credit (PTC) (\$/MW) for new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **ITCNew.csv**: Investment Tax Credit (ITC) (fraction) for new generation technologies. Columns are technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **ITCNewmy.csv**: Multi-year Investment Tax Credit (ITC) (fraction) for new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **ITCStorNewmy.csv**: Multi-year Investment Tax Credit (ITC) (fraction) for new storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **ITCOld.csv**: Investment Tax Credit (ITC) (fraction) for existing generation technologies. Columns are technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **ITCStorOld.csv**: Investment Tax Credit (ITC) (fraction) for existing storage technologies. Columns are technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

#### Bulk Resources (Utility and IPPGroup)

 - **DebtRatio.csv**: Debt ratio. Currently this file is only used by the IPPGroup and the debt ratio is defined over index\_p (e.g., one column called `ipp1`). Other supply-side agents (i.e., Utility, GreenDeveloper, DistributionUtility) have a debt ratio defined, but right now these are hard coded in those agents' constructors.

 - **COD.csv**: Cost of debt (%). Currently it is only used for IPP (for others they are constants). Column name is the IPP agent names in index\_p (e.g. `ipp1`).

 - **COE.csv**: Cost of equity (fraction). Currently it is only used for IPP (for others they are constants). Column name is the IPP agent names in index\_p (e.g. `ipp1`).

 - **Tax.csv**: Tax rate by IPP type. Column names are the IPP names in index\_p (e.g., `ipp1`).

 - **CumuTaxDepreOld.csv**: Cumulative tax depreciation of existing generation technologies. Columns are technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuTaxDepreOldmy.csv**: Multi-year cumulative tax depreciation of existing generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuTaxDepreStorOldmy.csv**: Multi-year cumulative tax depreciation of existing storage technologies, defined over years in index\_y. Columns are index\_y, followed by existing storage technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **CumuTaxDepreNew.csv**: Cumulative tax depreciation of new generation technologies. Columns are technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuTaxDepreNewmy.csv**: Multi-year cumulative tax depreciation of new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuTaxDepreStorNewmy.csv**: Multi-year cumulative tax depreciation of new storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **CumuAccoutDepreOld.csv**: Cumulative tax depreciation of existing generation technologies. Columns are technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuAccoutDepreNew.csv**: Cumulative accounting depreciation of new generation technologies. Columns are technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuAccoutDepreNewmy.csv**: Multi-year accounting depreciation of new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuAccoutDepreStorNewmy.csv**: Multi-year accounting depreciation of new storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **AnnualTaxDepreNewmy.csv**: Multi-year annual tax depreciation of new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **AnnualTaxDepreStorNewmy.csv**: Multi-year annual tax depreciation of new storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `gas-ct, gas-cc, upv, ...`).
 
 - **AnnualTaxDepreOldmy.csv**: Multi-year annual tax depreciation of existing generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **AnnualTaxDepreStorOldmy.csv**: Multi-year annual tax depreciation of existing storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **AnnualAccoutDepreNewmy.csv**: Multi-year annual accounting depreciation of existing generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **AnnualAccoutDepreStorNewmy.csv**: Multi-year annual accounting depreciation of new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **AnnualAccoutDepreOld.csv**: Annual accounting depreciation of existing generation technologies. Columns are technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuAccoutDepreOldmy.csv**: Multi-year accounting depreciation of existing generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuAccoutDepreStorOldmy.csv**: Multi-year accounting depreciation of existing storage technologies, defined over years in index\_y. Columns are index\_y, followed by storage technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **CumuITCAmortNew.csv**: Cumulative amortized ITC (fraction) for new generation technologies. Columns are technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuITCAmortNewmy.csv**: Multi-year cumulative amortized ITC (fraction) for new generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuITCAmortStorNewmy.csv**: Multi-year cumulative amortized ITC (fraction) for new storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **CumuITCAmortOldmy.csv**: Multi-year cumulative amortized ITC (fraction) for existing generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CumuITCAmortStorOldmy.csv**: Multi-year cumulative amortized ITC (fraction) for existing storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2 ...`).

 - **AnnualITCAmortNewmy.csv**: Multi-year annual amortized ITC (fraction) for new generation technologies, defined over years i index\_y. Columns are index\_y, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **AnnualITCAmortStorNewmy.csv**: Multi-year annual amortized ITC (fraction) for new storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **AnnualITCAmortOldmy.csv**: Multi-year annual amortized ITC (fraction) for existing generation technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **AnnualITCAmortStorOldmy.csv**: Multi-year annual amortized ITC (fraction) for existing storage technologies, defined over years in index\_y. Columns are index\_y, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).


#### Demand-Side (CustomerGroup)

 - **CustomerDebtRatio.csv**: Customer debt ratio by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (e.g., `index_z, Residential, Commercial, Industrial`).

 - **CustomerCOD.csv**: Customer cost of debt (fraction) by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (e.g., `index_z, Residential, Commercial, Industrial`).

 - **CustomerCOE.csv**: Customer cost of equity (fraction) by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (e.g., `index_z, Residential, Commercial, Industrial`).

 - **CustomerTax.csv**: Customer tax rate by customer type for zones in index\_z. Column names are index\_z, followed by customer group types in index\_h (e.g., `index_z, Residential, Commercial, Industrial`).

 - **DER_ITCNew.csv**: ITC for new DERs. Columns are BTM technology names in index\_m (e.g., `BTMPV, BTMStorage`).

 - **DER_ITCNewmy.csv**: Multi-year ITC for new DERs. Columns are index\_y, followed by BTM technology names in index\_m (e.g., `BTMPV, BTMStorage`).

#### SupplyChoice Use Case (GreenDeveloper)

*Agent code and data are nominally available, but have not been used in analysis, tested recently, or maintained.*

 - **RetailCompetition.csv**: Retail competition rates by year (used to calculate green tariff subscription percentages). Column names are years in index\_y (e.g., `2021, 2022`).

#### Rate-Making (Regulator)

 - **distribution_cost.csv**: Distribution cost (\$) for each zone in index\_z and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (e.g., `2021, 2022`). The filenames imply the corresponding type of cost.

 - **administration_cost.csv**: Administration cost (\$) for each zone in index\_z and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (e.g., `2021, 2022`). The filenames imply the corresponding type of cost.

 - **transmission_cost.csv**: Transmission cost (\$) for each zone in index\_z and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (e.g., `2021, 2022`). The filenames imply the corresponding type of cost.

 - **interconnection_cost.csv**: Interconnection cost (\$) for each zone in index\_z and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (e.g., `2021, 2022`). The filenames imply the corresponding type of cost.

 - **system_cost.csv**: System cost (\$ )for each zone in index\_z and for each year in index\_y. The column names are index\_z, followed by each year in index\_y (e.g., `2021, 2022`). The filenames imply the corresponding type of cost.

#### Bulk System Cost Parameters (Utility and IPPGroup)

 - **FOMNew.csv**: Fixed Operation and Maintenance cost (\$/MW-yr) of new generation technologies for each index\_z zone. Columns are index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FOMNewmy.csv**: Multi-year fixed Operation and Maintenance cost (\$/MW-yr) of new generation technologies for each zone in index\_z and each year in index\_y. Columns are index\_y, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FOMNewIPP.csv**: Fixed Operation and Maintenance cost (\$/MW-yr) of new generation technologies for each IPP in index\_p. Columns are index\_p, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FOMNewIPPmy.csv**: Multi-year fixed Operation and Maintenance cost (\$/MW-yr) of new generation technologies for each IPP in index\_p and year in index\_y. Columns are index\_y, index\_p, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapExNew.csv**: Capital cost (\$/MW) of new generation technologies for each index\_z zone. Columns are index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapExNewmy.csv**: Multi-year capital cost (\$/MW) of new generation technologies for each zone in index\_z and year in index\_y. Columns are index\_y, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapExNewIPP.csv**: Capital cost (\$/MW) of new generation technologies for each ipp in index\_p, for each zone in index\_z. Columns are index\_p, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapExNewIPPmy.csv**: Multi-year capital cost (\$/MW) of new generation technologies for each ipp in index\_p, zone in index\_z, and year in index\_y. Columns are index\_y, index\_p, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **StorCapExNew.csv**: Capital cost (\$/MW) of new storage technologies for each zone in index\_z. Columns are index\_z followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorCapExNewIPP.csv**: Capital cost (\$/MW) of new storage technologies for each IPP in index\_p, zone in index\_z. Columns are index\_p, index\_z followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorCapExNewIPPmy.csv**: Multi-year capital cost (\$/MW) of new storage technologies for each IPP in index\_p, zone in index\_z and year in index\_y. Columns are index\_y, index\_p,, index\_z followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapExOld.csv**: Capital cost (\$/MW) of existing generation technologies for each index\_z zone. Columns are index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapExStorOld.csv**: Capital cost (\$/MW) of existing storage technologies for each index\_z zone. Columns are index\_z, followed by storage technologies in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **FixedCostOld.csv**: Fixed cost (\$/MW) of existing generation technologies by zone in index\_z. Column names are index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FixedCostOldmy.csv**: Multi-year fixed cost (\$/MW ) of existing generation technologies by zone in index\_z, year in index\_y. Column names are index\_y, index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FixedCostOldIPP.csv**: Fixed cost (\$/MW) of existing generation technologies by IPP in index\_p and zone in index\_z. Column names are index\_p, index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FixedCostOldIPPmy.csv**: Multi-year fixed cost (\$/MW) of existing generation technologies by IPP in index\_p, zone in index\_z, year in index\_y. Column names are index\_y, index\_p, index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **FixedCostStorOld.csv**: Fixed cost (\$/MW ) of existing storage technologies by zone in index\_z. Column names are index\_z, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **FixedCostStorOldmy.csv**: Multi-year fixed cost (\$/MW ) of existing storage technologies by zone in index\_z, year in index\_y. Column names are index\_z, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **FixedCostStorOldIPP.csv**: Fixed cost (\$/MW) of existing storage technologies by IPP in index\_p and zone in index\_z. Column names are index\_p, index\_z, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **FixedCostStorOldIPPmy.csv**: Multi-year fixed cost (\$/MW) of existing storage technologies in by IPP in index\_p, zone in index\_z, and year in index\_y. Column names are index\_y, index\_p, index\_z, followed by technology names in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **VariableCostNew.csv**: Variable cost (\$/MW) of new generation technologies by zone in index\_z. Column names are index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostNewmy.csv**: Multi-year variable cost (\$/MW) of new generation technologies by zone in index\_z, year in index\_y. Column names are index\_y, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostOld.csv**: Variable cost (\$/MW) of existing generation technologies by zone in index\_z. Column names are index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostOldmy.csv**: Multi-year variable cost (\$/MW) of existing generation technologies by zone in index\_z, year in index\_y. Column names are index\_y, index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostNewIPP.csv**: Variable cost (\$/MW) of new generation technologies by IPP in index\_p, and zone in index\_z. Column names are index\_p, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostNewIPPmy.csv**: Multi-year variable cost (\$/MW) of new generation technologies by IPP in index\_p, zone in index\_z, and year in index\_y. Column names are index\_y, index\_p, index\_z, followed by technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostOldIPP.csv**: Variable cost (\$/MW) of existing generation technologies by IPP in index\_p and zone in index\_z. Column names are index\_p, index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **VariableCostOldIPPmy.csv**: Multi-year variable cost (\$/MW) of existing generation technologies by IPP in index\_p, zone in index\_z, and year in index\_y. Column names are index\_y, index\_p, index\_z, followed by technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **StorFOMNew.csv**: Fixed operation and maintenance cost (\$/MW-yr) of new storage technologies for each zone in index\_z. Columns are index\_z, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorFOMNewmy.csv**: Multi-year fixed operation and maintenance cost (\$/MW-yr) of new storage technologies for each zone in index\_z and year in index\_y. Columns are index\_y, index\_z, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorFOMNewIPP.csv**: Fixed operation and maintenance cost (\$/MW-yr) of new storage technologies for each IPP in index\_p and zone in index\_z. Columns are index\_p, index\_z, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorFOMNewIPPmy.csv**: Multi-year fixed operation and maintenance cost (\$/MW-yr) of new storage technologies for each IPP in index\_p, zone in index\_z and year in index\_y. Columns are index\_y, index\_p, index\_z, followed by technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

#### Customer Cost Parameters (CustomerGroup)

 - **CapExDER.csv**: Capital cost (\$/MW) of BTM technologies for each index\_z zone and customer type in index\_h. Columns are index\_z, index\_h, followed by technologies in index\_m (e.g., `BTMPV, BTMStorage`).

 - **CapExDERmy.csv**: Multi-year capital cost (\$/MW) of BTM technologies for each year in index\_y, zone in index\_z and customer type in index\_h. Columns are index\_y, index\_z, index\_h, followed by technologies in index\_m (e.g., `BTMPV, BTMStorage`).

 - **FOMDER.csv**: Fixed operation and maintenance cost (\$/MW-yr) of new BTM generation technologies for each zone in index\_z and for each customer type in index\_h. Columns are index\_z, index\_h, followed by technology names in index\_m (e.g., `BTMPV, BTMStorage`).

#### Capacity Market Parameters (IPPGroup)

 - **NetCONE.csv**: Net Cost of New Entry for capacity market (\$/MW-yr). Column names are the years in index\_y (e.g., `2021, 2022, ...`).

### Technology Parameters

#### Technology Performance for Existing and Future Resources (Utility and IPPGroup)

 - **StorRTENew.csv**: Round-trip efficiency (fraction) of new storage resources for each zone in index\_z. Column names are index\_z, followed by the technolgies in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorRTENewmy.csv**: Multi-year round-trip efficiency (fraction) of new storage resources by year in index\_y and zone in index\_z. Column names are index\_y, index\_z, followed by the technolgies in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorRTENewIPP.csv**: Multi-year round-trip efficiency (fraction) of new storage resources by year in index\_y, and zone in index\_z. Column names are index\_y, index\_z, followed by the technolgies in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorRTENewIPPmy.csv**: Multi-year round-trip efficiency (fraction) of new storage resources by year in index\_y, IPP in index\_p and zone in index\_z. Column names are index\_y, index\_p, index\_z, followed by the technolgies in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **StorRTEOldmy.csv**: Round-trip efficiency (fraction) of existing storage resources by zone in index\_z. Column names are index\_z, followed by the technolgies in index\_stor\_existing (e.g., `battery_1, battery_2, ...` etc.).

 - **StorRTEOldIPPmy.csv**: Multi-year round-trip efficiency (fraction) of exisiting storage resources by year in index\_y, IPP in index\_p and zone in index\_z. Column names are index\_y, index\_p, index\_z, followed by the technolgies in index\_stor\_existing (e.g., `battery_1, battery_2, ...`).

 - **NewStorDuration.csv**: Duration in hours of new storage technologies. Column names are the new storage technology names in index\_stor\_new (e.g., `battery_1, battery_2, ...`).

 - **NewStorInitialEnergy.csv**: Initial stored energy (Mwh) for new storage resources in index\_stor\_new by year in index\_y, zone in index\_z and day in index\_d. Column names are index\_y, index\_stor, index\_z, followed by the days in index\_d (e.g., `d001, d002, ...`).

 - **NewStorInitialEnergyIPP.csv**: Initial energy (MWh) for new storage resources in index\_stor\_new by year in index\_y, IPP in index\_p, zone in index\_z and day in index\_d. Column names are index\_y, index\_stor, index\_p, index\_z, followed by the days in index\_d (e.g., `d001, d002, ...`).

 - **CapacityCredit_new.csv**: Capacity credit (fraction) of new generation technologies by year in index\_y and zone in index\_z. Column names are index\_y, index\_z, followed by the existing technology names in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapacityCreditStor_new.csv**: Capacity credit (fraction) of new storage technologies by year in index\_y and zone in index\_z. Column names are index\_y, index\_z, followed by the new storage technology names in index\_stor_new (e.g., `battery_1, battery_2, ...`).

 - **CapacityCredit_old.csv**: Capacity credit (fraction) of existing generation technologies by year in index\_y and zone in index\_z. Column names are index\_y, index\_z, followed by the existing technology names in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **CapacityCreditStor_old.csv**: Capacity credit (fraction) of existing storage technologies by by year in index\_y and zone in index\_z. Column names are index\_y, index\_z, followed by the existing storage technology names in index\_stor\_existing  (e.g., `battery_1, battery_2, ...`).
 
 - **EmissionRateNewmy.csv**: Multi-year emission rate (pounds per million Btu) by year in index\_y, zone in index\_z and new technology types. Column names are index\_y, index\_z, followed by technolgies in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).

 - **EmissionRateNewIPPmy.csv**: Multi-year emission rate (pounds per million Btu) by year in index\_y, IPP in index\_p, zone in index\_z and new technology types. Column names are index\_y, index\_p, index\_z, followed by technolgies in index\_k\_new (e.g., `gas-ct, gas-cc, upv, ...`).
 
 - **EmissionRateOldmy.csv**: Multi-year emission rate (pounds per million Btu) by year in index\_y, zone in index\_z and new technology types. Column names are index\_y, index\_z, followed by technolgies in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

 - **EmissionRateOldIPPmy.csv**: Multi-year emission rate (pounds per million Btu) by year in index\_y, IPP in index\_p, zone in index\_z and new technology types. Column names are index\_y, index\_p, index\_z, followed by technolgies in index\_k\_existing (e.g., `gas-ct, gas-cc, upv, ...`).

#### Technology Performance for Customers (CustomerGroup)

 - **rte_dist_stor.csv**: Distributed storage round trip efficiency (fraction) for a customer type in index\_h and a zone in index\_z. Column names are index\_z, followed by the customer type values in index\_h (e.g., `Commercial, Residential, Industrial`).

 - **duration_dist_stor.csv**: Distributed storage duration (hours) for a customer type in index\_h and a zone in index\_z. Column names are index\_z, followed by the customer type values in index\_h (e.g., `Commercial, Residential, Industrial`).

 - **OptimalDER.csv**: Optimal DER capacity (MWh) of the BTM technologies for each each zone in index\_z and customer group in index\_h. Column names are index\_z, index\_h, folowed by the technology names in index\_m (e.g., `BTMPV, BTMStorage`).

 - **OptimalDERmy.csv**: Multi-year optimal DER capacity (MWh) of the BTM technologies for each year in index\_y, zone in index\_z and customer group in index\_h. Column names are index\_y, index\_z, index\_h, folowed by the technology names in index\_m (e.g., `BTMPV, BTMStorage`).

 - **initial\_energy\_dist\_stor.csv**: Initial ditributed storage energy (MWh) for each zone in index\_z, customer group in index\_h and day in index\_d. Column names are index\_z, index\_h, followed by days in index\_d (e.g., `d001, d002, ...`).

#### Technology Adoption Parameters (CustomerGroup)

 - **Bass_P.csv**: Coefficient of innovation (p) values for the Bass diffusion model (used to construct the PV adoption model) for each cutomer type in index\_h, zone in index\_z and BTM technology combination. Column names are index\_z, index\_h, followed by the values in index\_m (e.g., `BTMPV, BTMStorage`).

 - **Bass_Q.csv**: Coefficient of imitation (q) values for the Bass diffusion model (used to construct the PV adoption model) for each cutomer type in index\_h, zone in index\_z and BTM technology combination. Column names are index\_z, index\_h, followed by the values in index\_m (e.g., `BTMPV, BTMStorage`).

 - **DC_length.csv**: Capacity market demand curve length for each year in index\_y. Column names are years in index\_y(e.g., `2021, 2022`).

 - **WTP.csv**: Wllingness to pay (\$/MW) for each year in index\_y. Column names are the years in index\_y (e.g., `2021, 2022`).
