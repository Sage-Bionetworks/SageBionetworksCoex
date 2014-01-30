# returns the content of gc() (space report) and proc.time() (elapsed time report)
# as a single vector
# cum is a cumulative data frame of results.  If not null the answer will be added
# as a new row and the entire frame returned.  Thus one may call cum<-spaceTime(cum)
# to build up a data frame of results
spaceTime<-function(label, cum=NULL) {
	space <-gc()
	time <- proc.time()
	spaceAndTime <- t(as.matrix(c(space[1,], space[2,], time)))
	rownames(spaceAndTime)<-label
	if (is.null(cum)) {
		return(spaceAndTime)
	} else {
		return(rbind(cum, spaceAndTime))
	}
}

spaceTimePlot<-function(spaceTimeStats) {
	elapsed<-spaceTimeStats[,"elapsed"]
	delta<-elapsed[2:length(elapsed)]-elapsed[1:length(elapsed)-1]
	opar <- par(mfrow=c(3,1))
	pie(delta, col=rainbow(length(delta)), mfrow=c(2,1), new=T)
	title(paste("Breakdown of ", secToHMS(elapsed[length(elapsed)]-elapsed[1])))
	
	mem<-spaceTimeStats[,c(2,4,6)]+spaceTimeStats[,c(8,10,12)]
	colnames(mem)<-c("used (mb)", "GC trigger (mb)", "max used (mb)")
	
	plotMemCol(mem, 1, isNew=F)
	plotMemCol(mem, 3, isNew=F)
}

secToHMS<-function(sec) {
	hours = trunc(sec/3600)
	remainder = sec-hours*3600
	minutes = trunc(remainder/60)
	seconds = round(remainder-minutes*60)
	paste(paste(hours, "h", sep=""), paste(minutes, "m", sep=""), paste(seconds, "s", sep=""), sep=":")
}


plotMemCol<-function(mem, col, isNew) {
	# plot 'used (mb)'
	plot(1:(dim(mem)[1]), mem[,col], "s", xlab="step", ylab="mb", new=isNew) # draw lines
	text(1:(dim(mem)[1]), mem[,col], rownames(mem), cex=.8, pos=4, srt=-30)#label
	title(colnames(mem)[col])
}
