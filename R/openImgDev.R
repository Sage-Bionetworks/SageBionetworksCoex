openImgDev <- 
		function(imgname, iwidth = 1024, iheight = 1024, ipointsize = 12)
{

	## if (!require(lattice)) stop("Package lattice is required.")
	
	imgtype = .getFileExtension(imgname)
	
	if (imgtype=="ps"){
		postscript(imgname,width=iwidth, height=iheight, pointsize=ipointsize)
	}else if (imgtype=="png"){
		png(imgname, width=iwidth, height=iheight, pointsize=ipointsize)
	}else if (imgtype=="jpg" || imgtype=="jpeg"){
		jpeg(imgname, width=iwidth, height=iheight, pointsize=ipointsize,quality =100)
	}else if (imgtype=="pdf" || imgtype=="PDF"){
		pdf(imgname)
		return   
	}else{
		png(imgname, width=iwidth, height=iheight, pointsize=ipointsize)
	}
	trellis.device(new = FALSE, col = TRUE) 
}


.getFileExtension <- 
		function(fullfname)
{
	splitted=unlist( strsplit(fullfname, "\\.") )
	
	if( length(splitted) >1){
		return (splitted[length(splitted)])
	} else{
		return ("")
	}
}