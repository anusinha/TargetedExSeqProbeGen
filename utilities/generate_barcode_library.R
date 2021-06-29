library("DNABarcodes")

bcLength <- 7
minHammingDist <- 3
popSize <- 1000
numIter <- 500

mySet <- create.dnabarcodes(bcLength, dist = minHammingDist, heuristic = "ashlock", cores = 4, filter.triplets = TRUE, filter.gc = FALSE, population = popSize, iterations = numIter)
mySet <- paste0("'", gsub("T", "3", gsub("G", "2", gsub("C", "1", gsub("A", "0", mySet)))))
write(mySet, file = "barcode_set.csv")
