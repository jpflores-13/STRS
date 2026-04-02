# aggregateTAD.R --------------------------------------------------------- ## Aggregate Hi-C TADs from .hic files with scaling normalization

# Libraries ----
library(strawr)
library(data.table)
library(reshape2)
library(raster)

#' Aggregate TAD matrices from a .hic file
#'
#' Removes interchromosomal entries, buffers each TAD by a fraction of its
#' size, resamples each matrix to a common size, and normalizes by total
#' counts and number of TADs.
#'
#' @param loops data.frame. Regions to aggregate; first 6 columns must be in BEDPE format
#' @param hic character. Path to .hic file
#' @param buffer numeric. Fraction of TAD size to add on each side (default: 0.5)
#' @param res integer. Resolution of Hi-C data to use
#' @param size integer. Number of bins in aggregated matrix (default: 100)
#' @param norm character. Normalization method passed to strawr (default: "SCALE")
#' @return matrix. Normalized aggregated TAD matrix (dimensions: size x size)
aggregateTAD <- function(loops, hic, buffer = .5, res, size = 100, norm = "SCALE") {
  # remove interchrom loops
  loops = loops[which(loops[,1] == loops[,4]),]

  ## can modify to take Ginteractions

  # define window to plot
  loops$size = (loops[,5]- loops[,2])/res
  loops$bufferstart = res*(loops[,2]/res - round(loops$size*buffer))
  loops$bufferend = res*(loops[,5]/res + round(loops$size*buffer))
  loops$coords = paste0(loops[,1],":",loops$bufferstart,":",loops$bufferend)

  # create matrix to fill
  aggreTAD = matrix(0,nrow=size,ncol=size)

  # keep track of total counts
  totalcounts = 0

  # iterate through loops
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

    # check for bad loops
    if (min(rowSums(wideMat)) < 1 | min(colSums(wideMat < 1)))
    {
      print ("skipping")
      next()
    }

    # resize the matrix
    r_wideMat <- raster(wideMat)
    extent(r_wideMat) <- extent(c(-180, 180, -90, 90))

    resizedMat <- raster(ncol=size,  nrow=size)
    resizedMat <- resample(r_wideMat, resizedMat)

    # convert to matrix
    resizedMat = as.matrix(resizedMat)

    # update counts
    currentcounts = sum(resizedMat,na.rm = TRUE)
    totalcounts = totalcounts + currentcounts

    # normalize to more evenly weight different sized loops
    resizedMatNorm = resizedMat / currentcounts

    aggreTAD = aggreTAD + resizedMatNorm
  }

  # normalize for counts and number of loops
  aggreTAD = aggreTAD*totalcounts/nrow(loops)

  return (aggreTAD)
}

# Footer ----
message("aggregateTAD.R loaded. Available functions:")
message("  aggregateTAD(loops, hic, buffer, res, size, norm)")
message("Example: tad <- aggregateTAD(loops = mytads, hic = 'path/to/file.hic', res = 25000)")
