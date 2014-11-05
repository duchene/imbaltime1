beast2.path <- "~/PortableDucheneHo/BEAST2.1.1/bin/beast"
tree.ann.path <- "~/PortableDucheneHo/BEAST1.7.5/bin/treeannotator"
log.ann.path <- "~/PortableDucheneHo/BEAST1.7.5/bin/loganalyser"
code.path <- "~/PortableDucheneHo/apoint_scripts/"
pack.path <- "~/Desktop/sourceRlib"

require(phangorn)
require(geiger)
require(phytools)

dir.cur <- getwd()
setwd(code.path)
for(i in dir()) source(i)
setwd(dir.cur)
tree.sim <- read.tree(paste0('~/PortableDucheneHo/basetrees/het_', 'shallow','_','lowbal','.trees'))[[56]]
sim.dat <- molrate.sim.tree(tree = tree.sim, func = 'uncorlog')
sim.cal <- getCals(sim.dat, scheme = 'shallow', N = 4)
system(paste("mkdir", "uncorlog_lowbal_shallow_1"))
setwd("uncorlog_lowbal_shallow_1")
get.xml2het(tree = tree.sim, cal.data = sim.cal, sim.data = sim.dat, file.name = "uncorlog_lowbal_shallow_1")
beast.command <- paste("beast" ,paste0("uncorlog_lowbal_shallow_1", ".xml"))
## if beagle fails delete the -beagle command
system(beast.command)
write.tree(sim.dat[[1]], file = paste0("uncorlog_lowbal_shallow_1", "_sim.tre"))
write.table(sim.dat[[3]], file = paste0("uncorlog_lowbal_shallow_1", "_data_sim", ".csv"), sep = ",", row.names = F)
write.table(get.rate.descendant.pairs(sim.dat[[3]]), file = paste0("uncorlog_lowbal_shallow_1", "sim_rate_pairs.csv"), sep = ",", row.names = F)
dput(sim.dat, file = paste0("uncorlog_lowbal_shallow_1", ".Rdata"))
# This is redundant. Remove if it occupies too much memmory
treeannotator.command <- paste(tree.ann.path, "simulationduchene.trees", "estimated.tree")
system(treeannotator.command)
if(file.exists("estimated.tree") == F){
treesfile <- readLines("simulationduchene.trees")
open <- sapply(gregexpr("[(]", treesfile), length)
close <- sapply(gregexpr("[)]", treesfile), length)
ntrees <- sapply(gregexpr("tree", treesfile), length)
lastbit <- grepl("[;]$", treesfile)
treesfile <- treesfile[which(ntrees == 1 & open == 49 & close == 49 & lastbit)]
treesfile <- c("#NEXUS

Begin taxa;
	Dimensions ntax=50;
		Taxlabels
			s1 
			s2 
			s3 
			s4 
			s5 
			s6 
			s7 
			s8 
			s9 
			s10 
			s11 
			s12 
			s13 
			s14 
			s15 
			s16 
			s17 
			s18 
			s19 
			s20 
			s21 
			s22 
			s23 
			s24 
			s25 
			s26 
			s27 
			s28 
			s29 
			s30 
			s31 
			s32 
			s33 
			s34 
			s35 
			s36 
			s37 
			s38 
			s39 
			s40 
			s41 
			s42 
			s43 
			s44 
			s45 
			s46 
			s47 
			s48 
			s49 
			s50
			;
End;
Begin trees;
	Translate
		   1 s1,
		   2 s2,
		   3 s3,
		   4 s4,
		   5 s5,
		   6 s6,
		   7 s7,
		   8 s8,
		   9 s9,
		  10 s10,
		  11 s11,
		  12 s12,
		  13 s13,
		  14 s14,
		  15 s15,
		  16 s16,
		  17 s17,
		  18 s18,
		  19 s19,
		  20 s20,
		  21 s21,
		  22 s22,
		  23 s23,
		  24 s24,
		  25 s25,
		  26 s26,
		  27 s27,
		  28 s28,
		  29 s29,
		  30 s30,
		  31 s31,
		  32 s32,
		  33 s33,
		  34 s34,
		  35 s35,
		  36 s36,
		  37 s37,
		  38 s38,
		  39 s39,
		  40 s40,
		  41 s41,
		  42 s42,
		  43 s43,
		  44 s44,
		  45 s45,
		  46 s46,
		  47 s47,
		  48 s48,
		  49 s49,
		  50 s50
;
", treesfile, "End;")
writeLines(treesfile, con = "simulationduchene1.trees")
treeannotator.command <- paste(tree.ann.path, "simulationduchene1.trees", "estimated.tree")
system(treeannotator.command)
}
loganalyser.command <- paste(log.ann.path, paste0("uncorlog_lowbal_shallow_1", ".log"), paste0("uncorlog_lowbal_shallow_1", "_log.csv"))
system(loganalyser.command)
erase.files <- grep("[.]log|[.]trees", dir(), value=T)
sapply(erase.files, function(x) system(paste("rm", x)))
fol <- vector()
    setwd("..")
    for(i in length(dir())){
    if(file.exists(paste0(dir()[i], "/estimated.tree"))){ fol <- append(fol, dir()[i]) }}
    
    if(length(fol) >= 
135
){
    if(as.numeric(strsplit(getwd(), "hetrun")[[1]][2]) < 10){
    setwd(paste0("../hetrun", as.numeric(strsplit(getwd(), "hetrun")[[1]][2]) + 1))
    system("qsub submit_all_job.sh")
    }}
    
