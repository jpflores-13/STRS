#' @param loops data frame of regions to aggregate. Fist 6 columns must be in BEPDE format
#' @param hic  path to .hic file to use for aggregation
#' @param buffer pixels on either side
#' @param res resolution of hic data to use
#' @param size integer descibing the number of bins in the aggregated matix. (e.g. a size of 100 would result in an aggregated matrix with dimensions 100 x 100)

### Define a function that builds aggregate TADS
aggregateLoops <- function(loops, hic, buffer=.5, res, size=10, norm= "NONE")
{
  # remove interchrom loops
  loops = loops[which(loops[,1] == loops[,4]),]
  head(loops)
  # define window to plot
  loops$bufferstart = res*(loops[,2]/res - buffer)
  loops$bufferend = res*(loops[,5]/res + buffer)
  loops$coords = paste0(loops[,1],":",loops$bufferstart,":",loops$bufferend)
  
  # create matrix to fill
  aggreLoop = matrix(0,nrow=2*buffer+1,ncol=2*buffer+1)
  
  # keep track of total counts
  totalcounts = 0
  
  # iterate through loops
  successfulloops = 0
  for (i in 1:nrow(loops))
  {
    # print update
    print (paste(i, "of",nrow(loops),"loops"))
    
    # get loop info
    loop = loops[i,]
    
    # get pixels
    sparseMat = as.data.table(strawr::straw(norm, hic, loop$coords, loop$coords, "BP", res))
    
    # define bins
    startcoord = loop$bufferstart
    endcoord   = loop$bufferend
    bins <- seq(from = startcoord, to = endcoord, by = res)
    
    # make empty long format matrix
    longmat = as.data.table(expand.grid(bins,bins))
    longmat$counts = 0
    colnames(longmat) = c("x","y","counts")
    
    ## Set keys
    setkeyv(sparseMat, c('x', 'y'))
    
    ## Get counts by key
    longmat$counts <- sparseMat[longmat]$counts
    
    ## Set unmatched counts (NA) to 0
    longmat[is.na(counts), counts := 0]
    
    # convert to wide matrix
    wideMat <- reshape2::acast(longmat, x ~ y, value.var = 'counts')
    
    # make symmetric
    wideMat[lower.tri(wideMat)] = t(wideMat)[lower.tri(wideMat)]
    
    #subset to upper left
    wideMat = wideMat [1:(2*buffer+1),(dim(wideMat)[1] - (2*buffer):0)]
    
    # check for bad loops
    if (min(rowSums(wideMat)) < 1 | min(colSums(wideMat)) < 1 )
    {
      print (min(colSums(wideMat)))
      print (min(rowSums(wideMat)))
      print ("skipping")
      next()
    }
    successfulloops = successfulloops + 1
    aggreLoop = aggreLoop + wideMat
  }
  print (paste(successfulloops, "good loops"))
  return (aggreLoop/successfulloops)
}