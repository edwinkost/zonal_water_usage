

table_filename <- 
table = read.table(table_filename, header=FALSE, colClasses = "character", comment.char = "")

reg_linear <- lm( as.numeric(table[,2]) ~ as.numeric(table[,1]) - 1)
summary(reg_linear)



