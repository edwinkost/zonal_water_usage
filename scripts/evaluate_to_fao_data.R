
require(ggplot2)
require(foreign)

# FAO data
table_folder  = "/scratch-shared/edwinhs/country_water_use_for_pcrglobwb2.0_paper/table/"
variable_name = ""

Total water withdrawal
Fresh surface water withdrawal (primary and secondary)
Fresh groundwater withdrawal (primary and secondary)


# summary of PCR-GLOBWB output
pcrglobwb_output_folder = "/scratch-shared/edwinhs/country_water_use_for_pcrglobwb2.0_paper/table_summary/"
dir.create(output_folder)

# evaluation/validation will be done for the following periods
periods_in_year = seq(1960, 2015, 5)

for (i_period in 1:length(periods_in_year)){

mid_year = periods_in_year[i_period]
sta_year = mid_year - 2
end_year = mid_year + 2

# PCR-GLOBWB runs until 2015 only
if (mid_year == 2015) {end_year = 2015}

# loop over FAO country id:
# - load FAO data
# - find the corresponding WB id and load PCR-GLOBWB data

# save it on table

} # end for loop for i_period

# plot

2015 2010 2005 
2000 1995 1990 
1985 1980 1975 
1970 1965 1960

# plots: two figures: 
# - total, surface water, groundwater 
# - irrigation + livestock, industry, domestic
# - 


# merge the aforementioned variables to the final table
final_table = cbind(selected_table, 
                    average, std_dev, cor_pearson, cor_spearman, cor_kendall, cor_pearson__p_value, cor_spearman_p_value, cor_kendall__p_value, 
                    lm_slope, lm_r_squared, lm_adj_r_squared,
                    avg_1958_1962_1960, std_1958_1962_1960,
                    avg_1963_1967_1965, std_1963_1967_1965,
                    avg_1968_1972_1970, std_1968_1972_1970,
                    avg_1973_1977_1975, std_1973_1977_1975,
                    avg_1978_1982_1980, std_1978_1982_1980,
                    avg_1983_1987_1985, std_1983_1987_1985,
                    avg_1988_1992_1990, std_1988_1992_1990,
                    avg_1993_1997_1995, std_1993_1997_1995,
                    avg_1998_2002_2000, std_1998_2002_2000,
                    avg_2003_2007_2005, std_2003_2007_2005,
                    avg_2008_2012_2010, std_2008_2012_2010,
                    avg_2013_2015_2015, std_2013_2015_2015
                    )


# sort table 
final_table[order(as.numeric(as.character(final_table$FID))), ]

# write the final table to a txt file
output_file_name = paste(output_folder, "/", as.character(variable_names[i_variable]), ".txt", sep = "")
print(output_file_name)
write.table(final_table, output_file_name, sep = ";", row.names = FALSE)

} # end for loop for i_variable  
