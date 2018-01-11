
load("Yb_freq_table_pol.RDATA")
load("Yb_window_table.RDATA")

load("BD_freq_table_pol.RDATA")
load("BD_window_table.RDATA")

#Calculates D-statistic With frequency data
d.stat.func = function(left,right,freq.table.pol) {
	window = freq.table.pol[c(min(which(freq.table.pol[,2]>=left)):max(which(freq.table.pol[,2]<=right))),]
	if (length(window)==8) {d.stat = "NaN"}	
	if (length(window)!=8) {
		d.stat = sum(((1-window[,6])*window[,7]*window[,8])-(window[,6]*(1-window[,7])*window[,8]))/sum(((1-window[,6])*window[,7]*window[,8])+(window[,6]*(1-window[,7])*window[,8]))
	}	
	return(d.stat)
} 

#Same but with outgroup frequencies
d.stat.func = function(left,right,freq.table.pol) {
	#recover()
	window = freq.table.pol[c(min(which(freq.table.pol[,2]>=left)):max(which(freq.table.pol[,2]<=right))),]
	if (length(window)==8) {d.stat = "NaN"}
	if (length(window)!=8) {
		d.stat = sum(((1-window[,6])*window[,7]*window[,8]*window[,5])-(window[,6]*(1-window[,7])*window[,8]*window[,5]))/sum(((1-window[,6])*window[,7]*window[,8]*window[,5])+(window[,6]*(1-window[,7])*window[,8]*window[,5]))
	}	
	return(d.stat)
} 

#Same but with outgroup frequencies and no windows
d.stat.func = function(left,right,freq.table.pol) {
	#recover()
	window = freq.table.pol[left:right,]
	d.stat = sum(((1-window[,6])*window[,7]*window[,8]*window[,5])-(window[,6]*(1-window[,7])*window[,8]*window[,5]))/sum(((1-window[,6])*window[,7]*window[,8]*window[,5])+(window[,6]*(1-window[,7])*window[,8]*window[,5]))
	return(d.stat)
} 

#Uses d.stat.func to get  D-stats for all windows in window.table 
d.stat.vec=c(1:nrow(window.table))
for (i in c(1:nrow(window.table))) {
	print(i)
	d.stat.vec[i] = d.stat.func(window.table[i,1],window.table[i,2],freq.table.pol)
	d.table = cbind(window.table,d.stat.vec)	
}

save(d.table,file="Yb_d_table.RDATA")
save(d.table,file="BD_d_table.RDATA")
rm.vec = which(d.table[,3]=="NaN")
d.table = d.table[-rm.vec,]
window.table = window.table[-rm.vec,]
d.table = apply(d.table,2,as.numeric)
window.table = apply(window.table,2,as.numeric)
####################################################################
#Window

window.left = seq(min(raw.freqs[,2]),max(raw.freqs[,2]),5000)
window.left = window.left+c(0,c(1:(length(window.left)-1)))
window.right = window.left+5000
window.table = cbind(window.left,window.right)

matt.data = read.table("abbababa_C-P-M.txt")
matt.data = matt.data[which(matt.data[,1]=="chr1"),]
window.table=cbind(matt.data[which(matt.data[,1]=="chr1"),2],(matt.data[which(matt.data[,1]=="chr1"),2]+4999))

#chrom window.table
chr.vec = c(as.character(c(1:20)),"Z")
temp.chr="chr1"
for (i in 2:length(chr.vec)) {
	temp.chr = c(temp.chr,paste("chr",chr.vec[i],sep=""))
}
chr.vec=temp.chr
chr.vec=c(chr.vec,"Yb")
chr.vec=c(chr.vec,"BD")

start=c(1:length(chr.vec))
end=c(1:length(chr.vec))
for (i in 1:length(chr.vec)) {
	start[i] = min(which(scf.d[,1]==chr.vec[i]))
	end[i] = max(which(scf.d[,1]==chr.vec[i]))	
}
window.table=cbind(start,end)

#scaffold window.table
scaf = scaff.matrix[,1]
scf.vec = unique(scaff.matrix[,1])
start=c(1:length(scf.vec))
end=c(1:length(scf.vec))
for (i in 1:length(scf.vec)) {
	start[i] = min(which(scaf==scf.vec[i]))
	end[i] = max(which(scaf==scf.vec[i]))	
}
window.table=cbind(start,end)
scf.index = scf.vec[-c(which((window.table[,2]-window.table[,1])==0))]
window.table = window.table[-c(which((window.table[,2]-window.table[,1])==0)),]
#scf.index = scf.index[-which((window.table[,2]-window.table[,1]+1)<50)]
#window.table = window.table[-which((window.table[,2]-window.table[,1]+1)<50),]

#RAD window table
RAD.windows = c(0,0)
for (i in 1:length(scf.index)) {
	snp.pos = scaff.matrix[c(window.table[i,1]:window.table[i,2]),2]
	diff.pos = diff(snp.pos)
	breaks = which(diff.pos>1000)
	right = c(breaks,length(snp.pos))
	left = c(1,breaks+1)
	new.windows = cbind(left,right)
	snp.diffs = right-left
	filts = which(snp.diffs>20)
	first = min(which(scaff.matrix[,1]==scf.index[i]))-1
	new.windows = new.windows[filts,]+first
	RAD.windows = rbind(RAD.windows,new.windows)
}
RAD.windows = RAD.windows[-1,]
window.table = RAD.windows




diffs = as.numeric(RAD.windows[,2])-as.numeric(RAD.windows[,1])

plot(snp.pos)
points(breaks,snp.pos[breaks],col="red",lwd=3)


plot(positions[,2]-positions[,1],window.table[,2]-window.table[,1])
plot(temp.pos[,2]-temp.pos[,1],temp.wind[,2]-temp.wind[,1])
####################################################################

par(mfrow=c(2,1))
plot(c(1:(nrow(d.table)-2)),d.table[-c(919:920),2]-d.table[-c(919:920),1],type="l")
plot(c(1:nrow(d.table)),d.table[,3],type="l")

####################################################################
rm.vec=0
for (i in 1:4) {
	rm.vec = c(rm.vec,which(freq.table.pol[,(i+4)]=="NaN"))
}
rm.vec=rm.vec[-1]


png(file="~/kronforst_lab/figures/RAD_D.png",width=550,height=550)
plot(c(1:nrow(d.table)),d.table[,3],ylab="D",xlab="scaffold")
dev.off()






