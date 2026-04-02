# aggregateLoops.R ------------------------------------------------------- ## Aggregate Hi-C loops from .hic files

# Libraries ----
library(strawr)
library(data.table)
library(reshape2)

#' Aggregate loop matrices from a .hic file
#'
#' @param loops data.frame. Regions to aggregate; first 6 columns must be in BEDPE format
#' @param hic character. Path to .hic file
#' @param buffer numeric. Pixels on either side of loop anchor (default: 0.5)
#' @param res integer. Resolution of Hi-C data to use
#' @param size integer. Number of bins in aggregated matrix (default: 10)
#' @param norm character. Normalization method passed to strawr (default: "NONE")
#' @return matrix. Normalized aggregated loop matrix (dimensions: (2*buffer+1) x (2*buffer+1))
aggregateLoops <- function(loops, hic, buffer = .5, res, size = 10, norm = "NONE") {
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

# Footer ----
message("aggregateLoops.R loaded. Available functions:")
message("  aggregateLoops(loops, hic, buffer, res, size, norm)")
message("Example: agg <- aggregateLoops(loops = myloops, hic = 'path/to/file.hic', res = 10000)")
