args <- commandArgs(trailingOnly = TRUE)
	file <- c(args[1])
	headr <- c(args[2])
	if(headr=="T"){
			df <- read.csv(file, header = T, check.names = F)
			newdf <- aggregate(. ~ Species, df, FUN = sum)

	}else{
			df <- read.csv(file, header = F)
			colnames(df) <- c("Species","Abundance")
			newdf <- aggregate(. ~ Species, df, FUN = sum)
	}
	write.table(newdf,file,row.names = F,quote = F)
