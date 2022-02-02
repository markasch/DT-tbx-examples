library("pcalg") 
data("gmG") # load pre-defined data
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x)) 
pc.gmG <- pc(suffStat, indepTest = gaussCItest, p = ncol(gmG8$x), alpha = 0.01) 
stopifnot(require(Rgraphviz))# needed for all our graph plots 
par(mfrow = c(1,2)) # plot true and estimated DAGs
plot(gmG8$g, main = "") ; plot(pc.gmG, main = "")
## Data generation
## - used to generate "gmG"     
set.seed(40)     
p <- 8     
n <- 5000     
## true DAG:     
vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))     
gGtrue <- randomDAG(p, prob = 0.3, V = vars)     
gmG  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)     
gmG8 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)