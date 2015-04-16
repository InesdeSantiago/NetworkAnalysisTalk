############################################
#Example analysis with the RTN package
############################################

library(RTN)
data(dt4rtn)

#Compute network
rtni <- new("TNI", gexp=dt4rtn$gexp,transcriptionFactors=dt4rtn$tfs)
rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)
rtni<-tni.dpi.filter(rtni)

#get results
tni.get(rtni,what="summary") #filtered network
refnet<-tni.get(rtni,what="refnet") #unfiltered network
tnet<-tni.get(rtni,what="tnet")

#plot number of targets per TF 
nr_targets <- apply(tnet, 2, function (x) sum(x != 0))
summary(nr_targets)
barplot(sort(nr_targets), las=2, cex.names=0.5, 
xlab="TFs (probe ids)", ylab="Number of targets", cex.lab=3,
main="Number of targets per TF (DPI network)", cex.main=3, cex.axis=2)

#MRA
rtna<-tni2tna.preprocess(object=rtni,
phenotype=dt4rtn$pheno,
hits=dt4rtn$hits,
phenoIDs=dt4rtn$phenoIDs
)

rtna<-tna.mra(rtna) #hypergeometric test
rtna<-tna.gsea1(rtna) #GSEA test

#get results
tna.get(rtna,what="mra")
tna.get(rtna,what="gsea1")

save(rtni, rtna, file="/Users/santia01/Desktop/example_data.Rda")

############################################
#RedeR plot
############################################

getJC <- function(ref,minRegulonSize=20)
	{
	library(arules)
	ref <- ref[,colSums(ref != 0)>=minRegulonSize]
	ref[ref != 0] <- 1 #Remove weights from edges
	d_jaccard <- dissimilarity(t(ref))
	1 - as.matrix(d_jaccard)
	}

library(RedeR)
library(RColorBrewer)	
refnet<-tni.get(rtni,what="refnet") #unfiltered network
cormat <- getJC(refnet)

rdp<-RedPort()
calld(rdp,maxlag=5000)
resetd(rdp)
g=graph.adjacency(cormat,mode="undirected",weighted=TRUE,diag = FALSE)
relax(rdp, 300, 50, 100, ps = TRUE)
addGraph(rdp, g, layout = NULL)

#Add node symbols
tfs<-tni.get(rtni, what = "tfs")
idx<-tfs%in%colnames(cormat)
masterids<-data.frame(probeid=tfs[idx],symbol=names(tfs[idx]))
rownames(masterids)<-tfs[idx]

#get regulon's degree
bin <-tni.get(rtni,what="refnet")
bin[bin != 0]=1
degree<-colSums(bin)
masterids<-data.frame(masterids, degree=degree[rownames(masterids)]) 

#add color to nodes that are enriched (GSEA test)
regs <- tna.get(rtna,what="gsea1")$Regulon
hits <- masterids$symbol %in% regs
masterids<-data.frame(masterids, HitsGen = as.numeric(hits))

# set attributes  
g<-att.mapv(g=g,dat=masterids,refcol=1)
g<-att.setv(g=g, from="symbol", to='nodeAlias')
g<-att.sete(g=g, from="weight", to='edgeWidth', nquant=10, xlim=c(5,100,1))
g<-att.setv(g=g, from="degree", to='nodeSize', xlim=c(50,210,1), nquant=10, roundleg=1)
V(g)$nodeFontColor="black"
E(g)$linkType="notnested"
mypalette<-colorRampPalette(brewer.pal(9,"Purples"), bias=1)(100)
E(g)$edgeColor<-mypalette[30]
V(g)$nodeLineColor<-mypalette[30]
g<-att.setv(g=g, from="HitsGen", to='nodeColor',cols=c("grey","orange"))
V(g)$nodeFontSize<-80

# add graph
addGraph(rdp, g, layout = NULL)
relax(rdp, 300, 50, 100, ps = TRUE)

# add graph
addLegend.color(rdp,g, position="topright", title="Enrichment", dxtitle=10, bend=0.6, ftsize=7, dxborder=5)
addLegend.size(rdp,g,position="topleft", vertical=TRUE, title="Regulon size", dxtitle=36, ftsize=7, dxborder=10, dyborder=10)
addLegend.size(rdp,g,position="topleft", type="edge", title="Jaccard coefficient", dxtitle=40, ftsize=7, dyborder=300, 
                 dxborder=10, intersp=40, edgelen=150)
      
relax(rdp, 300, 50, 100, ps = TRUE)
  