# Co_occurrence_Network.R
# -------------------------------------------------------------------------
# Construction of microbial co-occurrence networks based on Spearman rank
# correlation. Abundance profiles are filtered by minimum relative
# abundance and occurrence frequency. Pairwise correlations are computed
# with Benjamini-Hochberg correction; edges are retained where adjusted
# p < 0.05 and |rho| >= 0.7. The resulting adjacency matrix is converted
# to an undirected weighted network using the igraph package (Csardi &
# Nepusz, 2006) and annotated with taxonomic information.
#
# Input:
#   otu_abundance.csv    — OTU/ASV abundance table (rows = taxa,
#                          columns = samples, first column = row names)
#   taxonomy.csv         — taxonomic annotation table (rows = taxa,
#                          columns = Phylum, Class, Order, Family, Genus,
#                          Species; first column = row names matching the
#                          abundance table)
#
# Output:
#   otu_corr_matrix.txt       — filtered correlation matrix
#   network_adj_matrix.txt    — adjacency matrix of the final network
#   network_edge_list.txt     — edge list with correlation and direction
#   network_node_list.txt     — node list with taxonomic attributes
#   network.gml               — GML format for visualisation in Gephi
#
# Dependencies:
#   install.packages(c("Hmisc", "psych", "igraph"))
# -------------------------------------------------------------------------

library(Hmisc)
library(psych)
library(igraph)

# =========================================================================
# 1. Read input data
# =========================================================================
otu <- read.csv("otu_abundance.csv", row.names = 1, check.names = FALSE)
tax <- read.csv("taxonomy.csv",      row.names = 1, check.names = FALSE)

# =========================================================================
# 2. Filter low-abundance and low-frequency taxa
# =========================================================================
# Retain taxa with cumulative relative abundance >= 0.005
abundance_sum <- rowSums(otu)
otu_filtered  <- otu[abundance_sum >= 0.005, ]

# Retain taxa present in at least 5 samples
occurrence        <- otu_filtered
occurrence[occurrence > 0] <- 1
otu_filtered <- otu_filtered[rowSums(occurrence) >= 5, ]

# =========================================================================
# 3. Spearman correlation with multiple-testing correction
# =========================================================================
corr_result <- rcorr(t(otu_filtered), type = "spearman")

r_mat <- corr_result$r
p_mat <- corr_result$P

# Benjamini-Hochberg adjustment
p_mat <- p.adjust(p_mat, method = "BH")
dim(p_mat) <- dim(r_mat)

# Retain edges with adjusted p < 0.05 and |rho| >= 0.7
r_mat[p_mat > 0.05 | abs(r_mat) < 0.7] <- 0
diag(r_mat) <- 0

write.table(data.frame(r_mat, check.names = FALSE),
            "otu_corr_matrix.txt", col.names = NA, sep = "\t", quote = FALSE)

# =========================================================================
# 4. Construct undirected weighted network
# =========================================================================
g <- graph.adjacency(r_mat, weighted = TRUE, mode = "undirected")
g <- simplify(g)
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

# Store original correlation as a separate edge attribute; use absolute
# values as edge weights (weights must be non-negative)
E(g)$correlation <- E(g)$weight
E(g)$weight      <- abs(E(g)$weight)
E(g)$direction   <- ifelse(E(g)$correlation > 0, 1, -1)

# =========================================================================
# 5. Annotate nodes with taxonomy
# =========================================================================
tax_matched      <- tax[as.character(V(g)$name), ]
V(g)$Phylum      <- tax_matched$Phylum
V(g)$Class       <- tax_matched$Class
V(g)$Order       <- tax_matched$Order
V(g)$Family      <- tax_matched$Family
V(g)$Genus       <- tax_matched$Genus
V(g)$Species     <- tax_matched$Species

# =========================================================================
# 6. Export network files
# =========================================================================

# 6a. Adjacency matrix
adj_matrix <- as.matrix(get.adjacency(g, attr = "correlation"))
write.table(data.frame(adj_matrix, check.names = FALSE),
            "network_adj_matrix.txt", col.names = NA, sep = "\t", quote = FALSE)

# 6b. Edge list
edge_pairs <- data.frame(as_edgelist(g), stringsAsFactors = FALSE)
edge_list  <- data.frame(
  source      = edge_pairs[[1]],
  target      = edge_pairs[[2]],
  weight      = E(g)$weight,
  correlation = E(g)$correlation,
  direction   = E(g)$direction
)
write.table(edge_list, "network_edge_list.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# 6c. Node list
node_list <- data.frame(
  label   = V(g)$name,
  Phylum  = V(g)$Phylum,
  Class   = V(g)$Class,
  Order   = V(g)$Order,
  Family  = V(g)$Family,
  Genus   = V(g)$Genus,
  Species = V(g)$Species
)
write.table(node_list, "network_node_list.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# 6d. GML format for Gephi
write.graph(g, "network.gml", format = "gml")
