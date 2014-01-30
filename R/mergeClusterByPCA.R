mergeClusterByPCA <-
		function(mdendro, genecluster, pccluster, pcnames)
{
	
	colcode.reduced<-genecluster
	geneclusize    = as.list(table(colcode.reduced))
	finalgenecolor = as.character(genecluster)
	pcclunames = names(table(pccluster))
	
	# get ordered gene cluster assignment
	fdispcolor = .getDisplayColorSequence(as.character(genecluster[mdendro$order]))
	orderedgeneclus= fdispcolor[fdispcolor!=""]
	
	#we need include grey module, but use different position to differentiate it from the ordinary clusters
	# so that the recombination will not consider the greay module
	orderedgeneclus = c(orderedgeneclus, "grey")
	
	geneclupos     = c(1:(length(orderedgeneclus)-1), 99999)
	names(geneclupos) <- orderedgeneclus
	
	# get ordered cluster-sizes
	orderedGeneclusize  = NULL
	for (each in orderedgeneclus){
		orderedGeneclusize  = c(orderedGeneclusize, geneclusize[[each]])
	}
	names(orderedGeneclusize)  <- orderedgeneclus
	
	for (each in pcclunames ){
		if (each =="grey"){
			next
		}
		
		# PCs belongs to the current cluster 
		sel.pcs      = as.character(pccluster) == each
		if (sum(sel.pcs)<=1){
			next
		}
		sel.pcnames  = pcnames[sel.pcs]
		no.selpcs    = length(sel.pcnames) 
		
		#get the PCs' sequential positions in the original dendrogram
		sel.pcpos    = NULL
		for (eachmerged in sel.pcnames){
			sel.pcpos  = c(sel.pcpos, geneclupos[[eachmerged]])
		}
		sel.order = order(sel.pcpos)
		
		# get the neighboring segments of the PCs, then w
		#> order.selPCpos: 30 31 32 36 37 38 39 40
		#> shift.selPCpos: 29 30 31 32 36 37 38 39
		#>  diff.selPCpos:  1  1  1  4  1  1  1  1
		
		order.selPCnames= sel.pcnames[sel.order]
		order.selPCpos  = sel.pcpos[sel.order]
		shift.selPCpos  = c(order.selPCpos[1]-2, order.selPCpos[c(1:(no.selpcs-1))])
		diff.selPCpos   = order.selPCpos - shift.selPCpos
		bool.startpos   = (diff.selPCpos !=1)   
		startpos        = c(1:no.selpcs)[bool.startpos]
		endpos          = c(startpos[-1]-1, no.selpcs)
		nosegments      = length(endpos)
		
		#we choose the cluster with maximal size as the module assigment for this SEGMENT of clusters
		for (iseg in c(1:nosegments) ){
			seg = c(startpos[iseg]:endpos[iseg])
			seg.pcnames = order.selPCnames[seg]
			if (length(seg.pcnames)==1){
				next
			}
			mergedcolor= .getNameOfMaxElement(orderedGeneclusize, seg.pcnames)
			#cat("merged color=", mergedcolor, "\n")
			for (eachmerged in seg.pcnames){
				#cat("   ", eachmerged, "\n")
				finalgenecolor = ifelse(finalgenecolor==eachmerged, mergedcolor, finalgenecolor)
			}
		}
		
	}
	
	# now we need sort the module assignment with new names by using an existing function
	# to do so, we need convert the current colors into numerals which are the input of this 
	# existing function
	retcolors = .reassignModuleNames(finalgenecolor)
	
	retcolors
}

#find the middle of each cluster and label the middle position with the corrsponding color
.getDisplayColorSequence <-
		function(colordered)
{
	mylen = length(colordered)
	colordered2 = c(colordered[1], colordered[1:(mylen-1)] )
	colordiff   = (colordered != colordered2)
	colordiff[1] = TRUE
	colordiff[mylen] = TRUE
	#mydispcolor = ifelse(colordiff==TRUE, colordered, "")
	mydispcolor = rep("", mylen)
	mytrueseq = c(1:mylen)[colordiff]
	for (i in c(1:(length(mytrueseq)-1)) ){
		midi =  (mytrueseq[i] + mytrueseq[i+1])/2
		mydispcolor[midi] = colordered[midi]
	}
	fdispcolor = ifelse(mydispcolor=="grey", "", mydispcolor)
	fdispcolor
}

.getNameOfMaxElement <- 
		function(mytable, selectednames)
{
	maxsize=0
	maxname=""
	for (each in selectednames){
		if (mytable[[each]] >maxsize){
			maxsize = mytable[[each]]
			maxname = each
		}
	}
	return (maxname)
}

# now we need sort the module assignment with new names by using an existing function
# to do so, we need convert the current colors into numerals which are the input of this 
# existing function
.reassignModuleNames <- 
		function(mycolorvect, minmodulesize=0)
{
	fgenecolor = as.character(mycolorvect)
	ztable = table(fgenecolor)
	zfinal = rep(0, length(fgenecolor))
	iclu   = 1 
	for (each in names(ztable)){
		if (each=="grey")
			next
		if (ztable[[each]] < minmodulesize)
			next
		
		zfinal = ifelse(fgenecolor==each, iclu, zfinal)
		iclu=iclu+1
	}
	
	retcolors=.assignModuleColor(labelpred=zfinal, minsize1=1, anameallmodules=T, auseblackwhite=FALSE)
	retcolors
}

.assignModuleColor <- 
		function(labelpred, minsize1=50, anameallmodules=FALSE, auseblackwhite=FALSE, startlabel=0) 
{
	# here we define modules by using a height cut-off for the branches
	#labelpred= cutree2(hiercluster,h=heightcutoff)
	#cat(labelpred)
	
	#"0", grey module doesn't participate color assignment, directly assigned as "grey"
	labelpredNoZero = labelpred[ labelpred >0 ]
	sort1=-sort(-table(labelpredNoZero))
	sort1
	modulename= as.numeric(names(sort1))
	modulebranch= sort1 > minsize1
	no.modules=sum(modulebranch)
	
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
		#colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow")
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
	
	
	#"grey" means not in any module;
	colorhelp=rep("grey",length(labelpred))
	if ( no.modules==0){
		print("No mudule detected\n")
	} else{
		if ( no.modules > length(colorcode)  ){
			print( paste("Too many modules \n", as.character(no.modules)) )
		}
		
		if ( (anameallmodules==FALSE) | (no.modules <=length(colorcode)) ){
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
			blackwhite=c("red", "black")
			for(i in c((maxcolors+1):no.modules)){
				if(auseblackwhite==FALSE){
					icolor=paste("module", as.character(i), sep="")
				}else{#use balck white alternatively represent extra colors, for display only
					#here we use the ordered label to avoid put the same color for two neighboring clusters
					icolor=blackwhite[1+(as.integer(modulename[i])%%2) ]
				}
				extracolors=c(extracolors, icolor)
			}
			
			#combine the true-color code and the extra colorcode into a uniform colorcode for 
			#color assignment
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


