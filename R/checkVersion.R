

checkVersion<-function(minimumMajorVersion, minimumMinorVersion) 
{
	majorVersion<-as.numeric(sessionInfo()$R.version$major)
	minorVersion<-as.numeric(sessionInfo()$R.version$minor)
	if (majorVersion<minimumMajorVersion | (majorVersion==minimumMajorVersion & minorVersion<minimumMinorVersion)) 
		stop(paste("Insufficient version of R: major version must be at least", 
						minimumMajorVersion, "and minor version at least", minimumMinorVersion))
}