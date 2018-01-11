
load("scf_table.RDATA")
load("freq_matrix_SureSelect.RDATA")
load("Mapp_scf.RDATA")
load("Null_scf.RDATA")
load("scf_map.RDATA")


scf.map = data.frame(matrix(nrow=nrow(freq.matrix),ncol=9))
index=1
Null.scf = 0
Mapp.scf = 0
for (i in 1:nrow(scf.table)) {
	index.b = index-1
	scf.matches = freq.matrix[which(freq.matrix[,1]==scf.table[i,4]),]
	if (length(scf.matches)==0) {print(paste("NULL",as.character(i)));Null.scf=c(Null.scf,i);next}
	print(i)
	Mapp.scf = c(Mapp.scf,i)
	if (length(scf.matches)!=8) {
		end.row=nrow(scf.matches)
		new.pos=scf.matches[,2]
		alleles=scf.matches[,c(3:4)]		
		freqs=scf.matches[,c(5:8)]	
		new.index = nrow(scf.matches)
	}
	if (length(scf.matches)==8) {
		end.row=1;new.pos=scf.matches[2]
		alleles=scf.matches[3:4]
		freqs=scf.matches[5:8] 
		new.index = 1
	}
	scf.map[c((index):((index.b)+end.row)),1] = scf.table[i,1]
	scf.map[c((index):((index.b)+end.row)),2] = scf.table[i,4]
	scf.map[c((index):((index.b)+end.row)),3] = (as.numeric(scf.table[i,2])-1)+as.numeric(new.pos)
	scf.map[c((index):((index.b)+end.row)),c(4:5)] = alleles
	scf.map[c((index):((index.b)+end.row)),c(6:9)] = as.numeric(freqs)
	index = index+new.index
}

scf.d = scf.map

scf.d[which(scf.map[,1]=="chr1_unmapped"),1] = "chr1"
scf.d[which(scf.map[,1]=="chr2_unmapped"),1] = "chr2"
scf.d[which(scf.map[,1]=="chr3_unmapped"),1] = "chr3"
scf.d[which(scf.map[,1]=="chr6_unmapped"),1] = "chr6"
scf.d[which(scf.map[,1]=="chr7_unmapped"),1] = "chr7"
scf.d[which(scf.map[,1]=="chr8_unmapped"),1] = "chr8"
scf.d[which(scf.map[,1]=="chr11_unmapped"),1] = "chr11"
scf.d[which(scf.map[,1]=="chr12_unmapped"),1] = "chr12"
scf.d[which(scf.map[,1]=="chr13_unmapped"),1] = "chr13"
scf.d[which(scf.map[,1]=="chr15_unmapped"),1] = "chr15"
scf.d[which(scf.map[,1]=="chr16_unmapped"),1] = "chr16"
scf.d[which(scf.map[,1]=="chr17_unmapped"),1] = "chr17"
scf.d[which(scf.map[,1]=="chr18_unmapped"),1] = "chr18"
scf.d[which(scf.map[,1]=="chr19_unmapped"),1] = "chr19"
scf.d[which(scf.map[,1]=="chr20_unmapped"),1] = "chr20"
scf.d[which(scf.map[,1]=="chrZ_unmapped"),1] = "chrZ"

og.freqs=as.numeric(matrix[,1])
og.cons = round(og.freqs)
og.pol = og.freqs
og.pol[which(og.cons==0)] = 1 - og.freqs[which(og.cons==0)]


#polarizing 
no.chrom.d.pol=matrix
no.chrom.d.pol[,c(2:4)]= abs(matrix[,1]-matrix[,c(2:4)])
no.chrom.d.pol[,1] = og.pol	
filler = no.chrom.d.pol
for (i in 1:ncol(filler)) {
	filler[,i] = 0
}
no.chrom.d.pol = cbind(filler,no.chrom.d.pol)
no.chrom.d.pol[,2] = as.numeric(no.chrom.d[,2])





#window.tables for divergence
ol = window.table[22:23,]
bg = window.table[-c(22:23),]
window.list = list("ol"=ol,"bg"=bg)


#get rid of list
table=scf.d.pol
matrix = (matrix(nrow=nrow(table),ncol=4))
for (j in 1:4) {
	matrix[,j] = as.numeric(paste(unlist(table[,j+4])))
}
blank.matrix=(matrix(nrow=nrow(table),ncol=4))
matrix = cbind(blank.matrix,matrix)
matrix[,2] = as.numeric(paste(unlist(table[,2])))


scp -r joelsmith@beast.uchicago.edu:kronforst_lab/agl_project/scf_d.RDATA scf_d.RDATA

factor = c(1:21)
plot.table = cbind(factor,d.table[1:21,3])
pdf(file="agl_plot.pdf",width=8,height=6,pointsize = 10)
plot(plot.table,col = rgb(0,0,1,.7),lwd=10,ylim=c(-.2,.2),xlim=c(1,23))
points(cbind(c(22:23),d.table[22:23,3]),col = rgb(0,1,0,.7),lwd=10,ylim=c(-.2,.2),xlim=c(1,23))
abline(h=0,lwd=2)
dev.off()


png(file="~/kronforst_lab/figures/agl_D.png",width=550,height=550)
plot(c(1:279),d.table[,3],ylim=c(-.7,.5),xlim=c(0,300),ylab="D",xlab="scaffold")
points(c(278:279),d.table[278:279,3],col="red",lwd=3,ylim=c(-.7,.5),xlim=c(0,300),ylab="")
abline(h=0)
dev.off()

png(file="~/kronforst_lab/figures/agl_SNPs.png",width=550,height=550)
plot(c(1:279),d.table[,2]-d.table[,1],ylab="# of SNPs",xlab="scaffold",xlim=c(0,300),ylim=c(0,80000))
points(c(278:279),d.table[278:279,2]-d.table[278:279,1],lwd=3,col="red",xlim=c(0,300),ylim=c(0,80000))
dev.off()


