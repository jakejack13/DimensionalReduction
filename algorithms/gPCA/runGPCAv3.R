#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This is a script designed to be called from matlab to run generalized PCA analysis.
#
# The syntax is: 
#
# Rscript /whatever/runGPCA.R inputname.mat outputname.mat
#
# This goes with runGPCAv3.m
#
# Luke Sjulson, 2016-10-27


# test if there is at least one argument: if not, return an error
if (length(args)<3) {
    stop("This script requires three arguments: the directory, the input filename, and the output filename", call.=FALSE)
  } 

# args[1] = '/workspace2/forAdrien/data/BF1-1_CPP_2014-10-08'
# args[2] = 'tempGPCA.mat'
# args[3] = 'temptemp.mat'

# cd to working directory
setwd(args[1])

# # for testing only
# sink(args[3])
# cat("writing to test file\n")
# sink()

library("R.matlab")
library("generalizedPCA")

inputList = readMat(args[2], fixNames=FALSE)

y <- matrix(unlist(inputList$GPCAinput), ncol = unlist(inputList$ncol), byrow = FALSE)
GPCAoutput = generalizedPCA(y, k = unlist(inputList$k), M = unlist(inputList$M), family = unlist(inputList$family), quiet = FALSE, max_iters = unlist(inputList$max_iters), normalize = inputList$normalize>0)


# write to disk
fid <- file(args[3], "wb")
writeMat(fid, inputList=inputList, GPCAoutput=GPCAoutput, fixNames=FALSE)
close(fid)

