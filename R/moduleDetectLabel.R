##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: automatically detect & label modules
##
##
##   nameallmodules=FALSE: label modules with all possible colors
##                 =TRUE:  when # of modules exceeds length(colorcode), we use false color 
##                          names to label the reamining modules
##
##   useblackwhite=FALSE: label as normal
##                =TRUE:  label extra modules by black and white alternatively
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#The input is an hclust object.
moduleDetectLabel <- 
		function(hiercluster,heightcutoff=0.5,minsize1=20, nameallmodules=FALSE, useblackwhite=FALSE, startlabel=1) 
{
	
	# here we define modules by using a height cut-off for the branches
	labelpred= cutree2(hiercluster,h=heightcutoff)
	sort1    =-sort(-table(labelpred))
	#sort1
	modulename   = as.numeric(names(sort1))
	modulebranch = sort1 > minsize1
	no.modules   = sum(modulebranch)
	
	## if (useNumberAsLabel){
	##     # now make cluster label
	##     #
	##     colorcode=NULL
	##     for (i in c(startlabel:(startlabel+no.modules-1)) ){
	##         ipad = patchZeros(i)
	##         colorcode= c(colorcode, ipad)
	##     }
	## }else{
		# now we assume that there are fewer than 10 modules
		colorcode=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
				"pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
				"midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
				"gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
				"forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
				"thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
				"deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
				"brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
				"brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
				"brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
				"gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
				"gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
				"gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
				"gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
				"gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
				"gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
				"gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
				"gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
				"gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
				"gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")
	## }
	
	# “grey" means not in any module;
	colorhelp=rep("grey",length(labelpred))
	if ( no.modules==0){
		print("No mudule detected\n")
	}
	else{
		if ( no.modules > length(colorcode)  ){
			print( paste("Too many modules \n", as.character(no.modules)) )
		}
		
		if ( (nameallmodules==FALSE) || (no.modules <=length(colorcode)) ){
			labeledModules = min(no.modules, length(colorcode) )
			for (i in c(1:labeledModules)) {
				colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
			}
			## if(!useNumberAsLabel){
				colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
			## }
		}else{#nameallmodules==TRUE and no.modules >length(colorcode)
			maxcolors=length(colorcode)
			labeledModules = no.modules
			extracolors=NULL
			blackwhite=c("black", "white")
			for(i in c((maxcolors+1):no.modules)){
				if(useblackwhite==FALSE){
					icolor=paste("module", as.character(i), sep="")
				}else{#use balck white alternatively represent extra colors, for display only
					icolor=blackwhite[1+(i%%2) ]
				}
				extracolors=c(extracolors, icolor)
			}
			allcolorcode=c(colorcode, extracolors)
			for (i in c(1:labeledModules)) {
				colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
			}
			
			## if(!useNumberAsLabel){
				colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
			## }
		}
	}
	colorhelp
}

