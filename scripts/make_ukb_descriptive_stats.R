# Load the test_phenos.csv file and generate some statistics

library(data.table)

out.file <- "test_statistics.csv"
pheno.file <- "test_phenos.csv"
d <- fread(pheno.file, data.table=FALSE)

# Update ethnicity based on coding: https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001
d$EthnicBG <- ifelse(is.na(d$ethnicity) | grepl("^-", d$ethnicity), "Unknown", 
		ifelse(grepl("^1", d$ethnicity), "White", 
		  ifelse(grepl("^2", d$ethnicity), "Mixed",
		    ifelse(grepl("^3", d$ethnicity), "Asian or Asian British", 
		      ifelse(grepl("^4", d$ethnicity), "Black or Black British", 
			ifelse(grepl("^5", d$ethnicity), "Chinese", "Other"))))))

get.stats <- function(dat, stats) {
    out <- data.frame(Variable=stats, N=NA, Mean=NA, SD=NA, 
		      Median=NA, Min=NA, Max=NA, stringsAsFactors=FALSE)
    for (v in stats) {
	i <- which(stats==v)
	x <- na.omit(dat[,v])
        if (all(x %in% c(0,1))) {
	    out$N[i] <- sum(x)
	} else {
	    out$N[i] <- length(x)
	    out$Mean[i] <- mean(x)
	    out$SD[i] <- sd(x)
	    out$Median[i] <- median(x)
	    out$Min[i] <- min(x)
	    out$Max[i] <- max(x)
	}
    }
    return(out)
}

# variables to include in the table
stats <- c("sex", "bp_med", "age", "hdl", "tg", "ldl")

# levels of ethnicity group stratification
eth.levels <- names(sort(table(d$EthnicBG), decreasing=TRUE))

# each object in the list is a statistic table for one of the ethnicity groups
stats.list.by.group <- lapply(eth.levels, 
			      function(e) {
				# collect statistics for each lifestyle trait variable and all samples
				x <- d[d$EthnicBG == e,]
				stats.list <- lapply(list(x[!is.na(x$gPC1) & !is.na(x$stst) & x$stst==1,],
				                          x[!is.na(x$gPC1) & !is.na(x$stst) & x$stst==0,],
				                          x[!is.na(x$gPC1) & !is.na(x$ltst) & x$ltst==1,],
				                          x[!is.na(x$gPC1) & !is.na(x$ltst) & x$ltst==0,],
				                          x[!is.na(x$gPC1) & !is.na(x$stst) | !is.na(x$ltst),]),
				                     function(dat) get.stats(dat, stats))

				# label each table in the list before combining
				prefixes <- c("STST1.", "STST0.", "LTST1.", "LTST0.", "ALL.")
				for (i in 1:length(stats.list)) {
				    names(stats.list[[i]])[-1] <- paste0(prefixes[i], names(stats.list[[i]])[-1])
				}
				
				# combine tables in the style of the analysis plan descriptive statistics table
				stats.table <- Reduce(function(x,y) merge(x, y, by="Variable", sort=FALSE), stats.list)
				stats.table <- cbind(EthnicBG=rep(e, nrow(stats.table)), stats.table)
				return(stats.table)
		      })

# combine all the individual tables together
out.table <- do.call("rbind", stats.list.by.group)

# save the data
write.csv(out.table, out.file, row.names=F)


