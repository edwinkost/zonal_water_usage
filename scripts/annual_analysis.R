
require(ggplot2)
require(foreign)

# output folder
output_folder = "/scratch-shared/edwinhs/country_water_use_for_pcrglobwb2.0_paper/table_summary/"
output_folder = "/scratch-shared/edwinhs/country_water_use_and_demand_for_pcrglobwb2.0_paper/water_demand_table_summary/"
dir.create(output_folder)

# table folder
table_folder  = "/scratch-shared/edwinhs/country_water_use_for_pcrglobwb2.0_paper/table/"
table_folder  = "/scratch-shared/edwinhs/country_water_use_and_demand_for_pcrglobwb2.0_paper/table/"

# years used
starting_year = 1958 
year_used = seq(starting_year, 2015, 1)

# get the variable names
first_year_table = read.table(paste(table_folder, "/summary_", as.character(starting_year), "-12-31.txt", sep = ""), header = T)
# - variable names starting from the column 4
variable_names = names(first_year_table)[4:length(first_year_table)]

# get the class from the shape file
shp_dbf_table = read.dbf("/projects/0/dfguu/users/edwin/data/country_shp_from_tianyi/World_Polys_High.dbf")
FID = as.character(seq(0, length(shp_dbf_table$WB_NAMES) - 1, 1))
shp_dbf_table = cbind(FID, shp_dbf_table)
shp_dbf_table = shp_dbf_table[,!names(shp_dbf_table)=="WB_RULES"]

# arranging tables 
for (i_variable in 1:length(variable_names)){

# making an initial table 
selected_table = shp_dbf_table
# add a column contains the variable name (only to ease the reading)
selected_table = cbind(variable_names[i_variable], selected_table)
names(selected_table)[1] <- as.character(variable_names[i_variable])

for (i_year in 1:length(year_used)){

# open file
complete_table = read.table(paste(table_folder, "/summary_", as.character(year_used[i_year]), "-12-31.txt", sep = ""), header = T)

# read values for this variable and
cropped_table = data.frame(as.character(complete_table$class_id),  
                           complete_table[, which(names(complete_table) == variable_names[i_variable])])
names(cropped_table)[1] <- "FID"                            
names(cropped_table)[2] <- as.character(year_used[i_year])

# merge to the selected table
selected_table = merge(selected_table, cropped_table, by = "FID", all.x = TRUE)

} # end for loop for i_year


# calculate some statistics values for the entire period
# - starting and end columns
sta_col = which(names(selected_table) == as.character(starting_year))
end_col = which(names(selected_table) == as.character(year_used[length(year_used)]))
# - average and standard deviation values
average = apply(selected_table[,sta_col:end_col], 1, mean)
std_dev = apply(selected_table[,sta_col:end_col], 1, sd)
# - auto correlations (correlation to the years)
x = NA; x = year_used
y = NA; y = selected_table[,sta_col:end_col]; y[is.na(y)] <- 0
cor_pearson  = apply(y, 1, function(y) cor(x, y, method = "pearson" ))
cor_spearman = apply(y, 1, function(y) cor(x, y, method = "spearman"))
cor_kendall  = apply(y, 1, function(y) cor(x, y, method = "kendall" ))
cor_pearson__p_value =  apply(y, 1, function(y) cor.test(x, as.numeric(y), method = "pearson" )$p.value)
cor_spearman_p_value =  apply(y, 1, function(y) cor.test(x, as.numeric(y), method = "spearman")$p.value)
cor_kendall__p_value =  apply(y, 1, function(y) cor.test(x, as.numeric(y), method = "kendall" )$p.value)
# - regression analysis (to the years)
lm_slope          = apply(y, 1, function(y) lm(y~x)$coefficients[2])
lm_r_squared      = apply(y, 1, function(y) summary( lm(y~x) )$r.squared)
lm_adj_r_squared  = apply(y, 1, function(y) summary( lm(y~x) )$r.squared)
lm_slope[is.na(average)] = NA

# calculate average and standard deviation values for certain periods

periods_in_year = seq(1960, 2015, 5)

for (i_period in 1:length(periods_in_year)){

mid_year = periods_in_year[i_period]
sta_year = mid_year - 2
end_year = mid_year + 2

# PCR-GLOBWB runs until 2015 only
if (mid_year == 2015) {end_year = 2015}

# calculate some statistics values for a certain period
# - starting and end columns
sta_col = which(names(selected_table) == as.character(sta_year))
end_col = which(names(selected_table) == as.character(end_year))
# - average and standard deviation values
average_for_this_period = apply(selected_table[,sta_col:end_col], 1, mean)
std_dev_for_this_period = apply(selected_table[,sta_col:end_col], 1, sd)

# average and standard deviation values
assign(paste("avg", as.character(sta_year), as.character(end_year), as.character(mid_year), sep = "_"), average_for_this_period)
assign(paste("std", as.character(sta_year), as.character(end_year), as.character(mid_year), sep = "_"), std_dev_for_this_period)

} # end for loop for i_year

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
