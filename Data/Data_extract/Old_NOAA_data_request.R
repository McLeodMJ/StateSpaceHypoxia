https://www.nwfsc.noaa.gov/data/api/v1/source/trawl.individual_fact/selection.csv?
filters=year>2003,year<2016
&variables=operation_dim$operation_type,
haul_latitude_dim$latitude_degree_units,
haul_latitude_dim$latitude_minute_units,
haul_latitude_dim$latitude_second_units,
haul_longitude_dim$longitude_degree_units,
haul_longitude_dim$longitude_minute_units,
haul_longitude_dim$longitude_second_units,
operation_dim$operation_start_date_whid,
operation_dim$operation_end_date_whid,
date_dim$day_of_month,
date_dim$month_of_year,
date_dim$year,
field_identified_taxonomy_dim$grp_reg_depth_category,
best_available_taxonomy_dim$order_40,
scientific_name,
common_name,
partition,
length_cm,
length_observation_detail_dim$measurement_procurement,
length_observation_detail_dim$observation_type

# project=Groundfish Slope and Shelf Combination Survey- focal point but if data mentioned about in other projects, would like that data as well

https://www.nwfsc.noaa.gov/data/api/v1/source/trawl.individual_fact/selection.csv?
best_available_taxonomy_dim$genus_70=Sebastes
&variables=trawl_id,scientific_name,common_name,length_cm,length_type,year,year_stn_invalid,o2_at_gear_ml_per_l_der,
https://www.nwfsc.noaa.gov/data/api/v1/source/trawl.individual_fact/selection.csv?
common_name=lingcod
&variables=trawl_id,common_name,length_cm,length_type,year,year_stn_invalid,o2_at_gear_ml_per_l_der

https://www.nwfsc.noaa.gov/data/api/v1/source/trawl.operation_haul_fact/selection.csv?
  filters=year>2003,year<2016
&variables=operation_dim$operation_type,
latitude_dim$latitude_degree_units,
latitude_dim$latitude_minute_units,
latitude_dim$latitude_second_units,
longitude_dim$longitude_degree_units,
longitude_dim$longitude_minute_units,
longitude_dim$longitude_second_units,
operation_dim$operation_start_date_whid,
operation_dim$operation_end_date_whid,
date_dim$day_of_month,
date_dim$month_of_year,
date_dim$year,
operation_dim$operation_start_date_whid,
operation_dim$operation_end_date_whid,
depth_m,
o2_at_gear_ml_per_l_der

save
