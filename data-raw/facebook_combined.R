################################################################################
## Prepare and load Facebook data example --> edgelist & nnodes               ##
## (ids contains original node numbers)                                       ##
################################################################################

if (!file.exists("data-raw/facebook_combined.txt")) {
  download.file(
    "https://snap.stanford.edu/data/facebook_combined.txt.gz",
    "data-raw/facebook_combined.txt.gz"
  )
  system("gzip -d data-raw/facebook_combined.txt.gz")
}

raw_edgelist <- read.table("data-raw/facebook_combined.txt") # read-in data

## Number nodes continuously starting with 1 ###################################
## Note: all "missing" ids are otherwise treated as present, which leads to
##       nodes with no connections (maybe not what we want)
ids <- sort(unique(c(raw_edgelist[, 1], raw_edgelist[, 2])))
edgelist <- matrix(0, nrow = nrow(raw_edgelist), ncol = 2)
for (i in 1:nrow(raw_edgelist)) {
  edgelist[i, 1] <- which(raw_edgelist[i, 1] == ids)
  edgelist[i, 2] <- which(raw_edgelist[i, 2] == ids)
}
rm("i")
rm("raw_edgelist")

nnodes <- length(ids)

## Generate adjacency for full network
adjacency <- matrix(0, nrow = nnodes, ncol = nnodes)
# upper triangle
adjacency[edgelist] <- 1
# lower triangle
adjacency[edgelist[, c(2, 1)]] <- 1

# facebook <- adjacency # store adjacency only
library("network")
facebook <- as.network.matrix(adjacency, directed = FALSE)

save(facebook, file = "data/facebook.RData")
