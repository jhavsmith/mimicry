

# computes mean pairwise divergence per snp
divergence.func = function(raw.freqs,window.table) {
	#recover()
	pairwise.div.a.b = (raw.freqs[,6]*(1-raw.freqs[,7]))+((1-raw.freqs[,6])*raw.freqs[,7])
	pairwise.div.a.c = (raw.freqs[,6]*(1-raw.freqs[,8]))+((1-raw.freqs[,6])*raw.freqs[,8])
	pairwise.div.b.c = (raw.freqs[,7]*(1-raw.freqs[,8]))+((1-raw.freqs[,7])*raw.freqs[,8])
	snp.pos = raw.freqs[,2]	
	for (i in c(1:nrow(window.table))) {
		if (length(which(snp.pos<=window.table[i,2]))!=0) {		
			snps = c(min(which(snp.pos>=window.table[i,1])):max(which(snp.pos<=window.table[i,2])))
		}	
		cp = sum(pairwise.div.a.b[snps])/length(snps)
		cm = sum(pairwise.div.a.c[snps])/length(snps)
		pm = sum(pairwise.div.b.c[snps])/length(snps)
		divergence.table[i,] = c(cm,pm,cp)
	}
	return(divergence.table)
}


#per base-pair
divergence.func = function(raw.freqs,window.table) {
	#recover()
	pairwise.div.a.b = (raw.freqs[,6]*(1-raw.freqs[,7]))+((1-raw.freqs[,6])*raw.freqs[,7])
	pairwise.div.a.c = (raw.freqs[,6]*(1-raw.freqs[,8]))+((1-raw.freqs[,6])*raw.freqs[,8])
	pairwise.div.b.c = (raw.freqs[,7]*(1-raw.freqs[,8]))+((1-raw.freqs[,7])*raw.freqs[,8])
	snp.pos = raw.freqs[,2]	
	for (i in c(1:nrow(window.table))) {
		if (length(which(snp.pos<=window.table[i,2]))!=0) {		
			snps = c(min(which(snp.pos>=window.table[i,1])):max(which(snp.pos<=window.table[i,2])))
		}	
		cp = sum(pairwise.div.a.b[snps])/(window.table[i,2]-window.table[i,1])
		cm = sum(pairwise.div.a.c[snps])/(window.table[i,2]-window.table[i,1])
		pm = sum(pairwise.div.b.c[snps])/(window.table[i,2]-window.table[i,1])
		divergence.table[i,] = c(cm,pm,cp)
	}
	return(divergence.table)
}


#per base-pair. no windows
divergence.func = function(raw.freqs,window.table,positions) {
	#recover()
	pairwise.div.a.b = (raw.freqs[,6]*(1-raw.freqs[,7]))+((1-raw.freqs[,6])*raw.freqs[,7])
	pairwise.div.a.c = (raw.freqs[,6]*(1-raw.freqs[,8]))+((1-raw.freqs[,6])*raw.freqs[,8])
	pairwise.div.b.c = (raw.freqs[,7]*(1-raw.freqs[,8]))+((1-raw.freqs[,7])*raw.freqs[,8])
	for (i in c(1:nrow(window.table))) {
		cp = sum(pairwise.div.a.b[window.table[i,1]:window.table[i,2]])/(positions[i,2]-positions[i,1])
		cm = sum(pairwise.div.a.c[window.table[i,1]:window.table[i,2]])/(positions[i,2]-positions[i,1])
		pm = sum(pairwise.div.b.c[window.table[i,1]:window.table[i,2]])/(positions[i,2]-positions[i,1])
		divergence.table[i,] = c(cm,pm,cp)
	}
	return(divergence.table)
}

#fixed differences. per basepair. no windows
divergence.func = function(raw.freqs,window.table,positions) {
	#recover()
	pairwise.cp = raw.freqs[,6] - raw.freqs[,7]
	pairwise.cm = raw.freqs[,6] - raw.freqs[,8]
	pairwise.pm = raw.freqs[,7] - raw.freqs[,8]
	for (i in c(1:nrow(window.table))) {
		cp = (length(which(pairwise.cp[window.table[i,1]:window.table[i,2]]==-1))+length(which(pairwise.cp[window.table[i,1]:window.table[i,2]]==1)))/(positions[i,2]-positions[i,1])
		cm = (length(which(pairwise.cm[window.table[i,1]:window.table[i,2]]==-1))+length(which(pairwise.cm[window.table[i,1]:window.table[i,2]]==1)))/(positions[i,2]-positions[i,1])
		pm = (length(which(pairwise.pm[window.table[i,1]:window.table[i,2]]==-1))+length(which(pairwise.pm[window.table[i,1]:window.table[i,2]]==1)))/(positions[i,2]-positions[i,1])
		divergence.table[i,] = c(cm,pm,cp)
	}
	return(divergence.table)
}

#fixed differences. per basepair. 
divergence.func = function(raw.freqs,window.table) {
	#recover()
	pairwise.cp = raw.freqs[,6] - raw.freqs[,7]
	pairwise.cm = raw.freqs[,6] - raw.freqs[,8]
	pairwise.pm = raw.freqs[,7] - raw.freqs[,8]
	snp.pos = raw.freqs[,2]	
	for (i in c(1:nrow(window.table))) {
		if (length(which(snp.pos<=window.table[i,2]))!=0) {		
			snps = c(min(which(snp.pos>=window.table[i,1])):max(which(snp.pos<=window.table[i,2])))
		}	
		cp = (length(which(pairwise.cp[snps]==-1))+length(which(pairwise.cp[snps]==1)))/(window.table[i,2]-window.table[i,1])
		cm = (length(which(pairwise.cm[snps]==-1))+length(which(pairwise.cm[snps]==1)))/(window.table[i,2]-window.table[i,1])
		pm = (length(which(pairwise.pm[snps]==-1))+length(which(pairwise.pm[snps]==1)))/(window.table[i,2]-window.table[i,1])
		divergence.table[i,] = c(cm,pm,cp)
	}
	return(divergence.table)
}



#cutoff for absolute value of d
cutoff = quantile(abs(d.table[,3]),probs=.95)
ol.rows = c(which(abs(d.table[,3])>cutoff))
ol.table = d.table[ol.rows,]
ol = ol.table[which(d.table[ol.rows,3]>0),1:2]
bg = d.table[-ol.rows,1:2]
#plot(c(1:nrow(d.table)),d.table[,3],xlim=c(0,500),ylim=c(-1,1))
#points(ol.rows,d.table[ol.rows,3],pch=19,col="red",xlim=c(0,2000),ylim=c(-1,1))
window.list = list("ol"=ol,"bg"=bg)

########################################################################
#map positions to  scaf windows scf.d.pol is freq.table.pol
positions = c(0,0)
for (i in 1:nrow(window.table)) {
	print(i)
	new.positions = c(freq.table.pol[window.table[i,1],2],freq.table.pol[window.table[i,2],2])
	positions = rbind(positions,new.positions)
}
positions=positions[-1,]

###########################################################################

#Counts number of snps in each window
snp.count.func = function(raw.freqs,window.table) {
	#recover()
	snp.count.vec = c(1:nrow(window.table))	
	snp.pos = raw.freqs[,2]	
	for (i in c(1:nrow(window.table))) {
		snp.count.vec[i] = length(min(which(snp.pos>=window.table[i,1])):max(which(snp.pos<=window.table[i,2])))
	}
	return(snp.count.vec)
}


#Counts number of snps. No windows
snp.count.func = function(window.table) {
	snp.count.vec = c(1:nrow(window.table))	
	snp.count.vec = window.table[,2] - window.table[,1] + 1
	return(snp.count.vec)
}

all.snps = snp.count.func(raw.freqs,window.table)
#d.table=d.table[-which(all.snps<50),]
#window.table=window.table[-which(all.snps<50),]

#snp counts
for (j in 1:3) {
	snp.count.table = c(1:nrow(window.list[[j]]))
	snp.count.table = snp.count.func(window.list[[j]])
	if (j==1) {neg.snps = snp.count.table}
	if (j==2) {pos.snps = snp.count.table}
	if (j==3) {bg.snps = snp.count.table}	
}

#snp counts
for (j in 1:2) {
	snp.count.table = c(1:nrow(window.list[[j]]))
	snp.count.table = snp.count.func(raw.freqs,window.list[[j]])
	if (j==1) {ol.snps = snp.count.table}
	if (j==2) {bg.snps = snp.count.table}
}

#divergence tables
for (j in 1:2) {
	divergence.table = matrix(nrow=nrow(window.list[[j]]),ncol=3)
	divergence.table = divergence.func(raw.freqs,window.list[[j]],positions)
	if (j==1) {ol = divergence.table}
	if (j==2) {bg = divergence.table}
}

#divergence tables
for (j in 1:3) {
	divergence.table = matrix(nrow=nrow(window.list[[j]]),ncol=3)
	divergence.table = divergence.func(raw.freqs,window.list[[j]],positions)
	if (j==1) {neg.ol = divergence.table}
	if (j==2) {pos.ol = divergence.table}
	if (j==3) {bg = divergence.table}
}

########################################################################
#barplots
error.bar.func = function(data){
	#recover()
	sd.vec = apply(data,2,sd)	
	se.vec = sd.vec/sqrt(nrow(data))
	mean.vec = apply(data,2,mean)
	uppers = mean.vec+se.vec
	lowers = mean.vec-se.vec
	error.table = cbind(mean.vec,uppers)
	error.table = cbind(error.table,lowers)
	return(error.table)
}

#######################################################for outlier and bg d. 
ol.plot.table = error.bar.func(ol)
bg.plot.table = error.bar.func(bg)
all.means = ol.plot.table[,1]
all.means = rbind(all.means,bg.plot.table[,1])

bar.pos = barplot(all.means,beside=T)
x.vec = cbind(as.vector(bar.pos),as.vector(bar.pos))
temp.ups = ol.plot.table[,2]
temp.ups = rbind(temp.ups,bg.plot.table[,2])
ups = as.vector(temp.ups)
temp.lows = ol.plot.table[,3]
lows = as.vector(rbind(temp.lows,bg.plot.table[,3]))
bounds = cbind(ups,lows)

ups.t = ups
lows.t = lows
all.means.t = all.means
bounds.t = bounds
x.vec.t = x.vec

pdf(file="~/kronforst_lab/RAD_2.pdf",width=6,height=6,pointsize=10)
#par(mfrow=c(1,2))
colnames(all.means.t) = c("agl-tim","ama-tim","agl-ama") 
barplot(all.means.t,beside=T,ylim=c(0,.009),col=c("orange","darkolivegreen"),main="RAD,.2",ylab="divergence")
for (i in 1:nrow(x.vec.t)) {
	lines(x=x.vec.t[i,],y=bounds.t[i,])
	lines(x=c(x.vec.t[i,1]-.2,x.vec.t[i,1]+.2),y=cbind(ups.t,ups.t)[i,])
	lines(x=c(x.vec.t[i,1]-.2,x.vec.t[i,1]+.2),y=cbind(lows.t,lows.t)[i,])
}
#legend(6.5,.0065,c(".05 outliers","background"),c("orange","darkolivegreen"))
#plot(d.table[,3],ylim=c(-1,1),xlim=c(0,2000),lwd=1,ylab="D",xlab="locus",main="RAD")
#points(which(d.table[,3]>cutoff),d.table[which(d.table[,3]>cutoff),3],pch=20,lwd=4,col="red",ylim=c(-1,1),xlim=c(0,2000))	
#abline(h=0)
dev.off()

scp -r joelsmith@beast.uchicago.edu:kronforst_lab/RAD_posD.pdf RAD_posD.pdf

########################################################for pos and neg d
neg.ol.plot.table = error.bar.func(neg.ol)
pos.ol.plot.table = error.bar.func(pos.ol)
bg.plot.table = error.bar.func(bg)
all.means = rbind(pos.ol.plot.table[,1],neg.ol.plot.table[,1])
all.means = rbind(all.means,bg.plot.table[,1])

bar.pos = barplot(all.means,beside=T)
x.vec = cbind(as.vector(bar.pos),as.vector(bar.pos))
temp.ups = rbind(pos.ol.plot.table[,2],neg.ol.plot.table[,2])
temp.ups = rbind(temp.ups,bg.plot.table[,2])
ups = as.vector(temp.ups)
temp.lows = rbind(pos.ol.plot.table[,3],neg.ol.plot.table[,3])
lows = as.vector(rbind(temp.lows,bg.plot.table[,3]))
bounds = cbind(ups,lows)
	
ups.t = ups
lows.t = lows
all.means.t = all.means
bounds.t = bounds
x.vec.t = x.vec


pdf(file="~/kronforst_lab/RAD_95.pdf",width=12,height=6,pointsize=10)
op = par(mfrow=c(1,2))
colnames(all.means.t) = c("agl-tim","ama-tim","agl-ama") 
barplot(all.means.t,beside=T,ylim=c(0,.01),col=c("slate blue","orange","darkolivegreen"),main="RAD .609 cutoff",ylab="divergence")
legend(3,.01,c("pos d (ama <--> tim)","neg d (agl <--> tim)","background"),c("slate blue","orange","darkolivegreen"))	
for (i in 1:nrow(x.vec.t)) {
	lines(x=x.vec.t[i,],y=bounds.t[i,])
	lines(x=c(x.vec.t[i,1]-.2,x.vec.t[i,1]+.2),y=cbind(ups.t,ups.t)[i,])
	lines(x=c(x.vec.t[i,1]-.2,x.vec.t[i,1]+.2),y=cbind(lows.t,lows.t)[i,])
}
plot(c(1:nrow(d.table)),d.table[,3],xlim=c(0,2000),ylim=c(-1,1),main=".609 cutoff",ylab="D",xlab="locus")
points(ol.rows,d.table[ol.rows,3],pch=19,col="red",xlim=c(0,2000),ylim=c(-1,1))
abline(h=0)
dev.off()


load("RAD_freq_table_pol.RDATA")
load("RAD_window_table.RDATA")
load("RAD_d_table.RDATA")
raw.freqs=freq.table.pol


neg.div = neg.ol
pos.div = pos.ol
bg.div = bg









