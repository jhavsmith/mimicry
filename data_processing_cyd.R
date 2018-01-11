

# make a vector of chromosome labels
chr.vec = c(as.character(c(1:20)),"Z")
temp.chr="chr1"
for (i in 2:length(chr.vec)) {
	temp.chr = c(temp.chr,paste("chr",chr.vec[i],sep=""))
}
chr.vec=temp.chr
chrom.numbers = as.character(c(1:21))


#Filtering
hel.list = c("c511.a","c511.b","c512.a","c512.b","c513.a","c513.b","c514.a","c514.b","c515.a","c515.b","c563.a","c563.b","c614.a","c614.b","c630.a","c630.b","c639.a","c639.b","c640.a","c640.b","h665.a","h665.b","i02-210.a","i02-210.b","m523.a","m523.b","m524.a","m524.b","m525.a","m525.b","m589.a","m589.b","m675.a","m675.b","m676.a","m676.b","m682.a","m682.b","m683.a","m683.b","m687.a","m687.b","m689.a","m689.b","p516.a","p516.b","p517.a","p517.b","p518.a","p518.b","p519.a","p519.b","p520.a","p520.b","p591.a","p591.b","p596.a","p596.b","p690.a","p690.b","p694.a","p694.b","p696.a","p696.b")

for (j in 10:21) {
	data.list = get.data.func(j)
	data.list=biallele.func(data.list)
	data.list=missing.data.func(data.list)	
	save(data.list,file=paste("data_list",as.character(j),".RDATA",sep=''))
}
	
#Getting frequencies
for (j in 1:21) {
	load(paste("data_list",chrom.numbers[j],".RDATA",sep=''))
	freq.func(data.list)
}

#Getting sample sizes
for (j in 1:21) {
	print(chr.vec[j])
	load(paste("data_list",chrom.numbers[j],".RDATA",sep=''))
	n.func(data.list)
	print("finished n.func")
}

#Getting heterozygotes
for (j in 1:21) {
	print(chr.vec[j])
	load(paste("data_list",chrom.numbers[j],".RDATA",sep=''))
	het.func(data.list)
	print("finished het.func")
}


###############################################################################


get.data.func = function(j) {
	print(paste(as.character(j), "starting"))
	system(paste("cat 32-butterflies_a | grep ", chr.vec[j], " -w | cut -f 1,2,4,5,10-41 > ref_table",as.character(j),".txt",sep = ""))
	ref.table = read.table(paste("ref_table",as.character(j),".txt",sep=''))
	print("chr.table starting")
	temp.table = data.frame(matrix(nrow=length(ref.table[[1]]),ncol=4))
		temp.table[,2] = ref.table[,2]
		temp.table[,1] = unlist(lapply(ref.table[[1]],as.character))
		temp.table[,3] = unlist(lapply(ref.table[[3]],as.character))
 		temp.table[,4] = unlist(lapply(ref.table[[4]],as.character))
		chr.table = temp.table
	save(chr.table,file=paste("chr_table",as.character(j),".RDATA",sep=''))
	print("chr.table finished")
	print("gt.table starting")
	gt.list = ref.table[,c(5:36)]
	i=1
	gt.temp = unlist(gt.list[,i],toString)
	gt.temp = paste(gt.temp,collapse=" ")
	gt.temp = strsplit(gt.temp," ")
	gt.vec = gt.temp[[1]]
	allele.a = substring(gt.vec,1,1)
	allele.b = substring(gt.vec,3,3)
	gt.table = cbind(allele.a,allele.b)
	for (i in c(2:32)) {
		print(i)
		gt.temp = unlist(gt.list[,i],toString)
		gt.temp = paste(gt.temp,collapse=" ")
		gt.temp = strsplit(gt.temp," ")
		gt.vec = gt.temp[[1]]
		allele.a = substring(gt.vec,1,1)
		allele.b = substring(gt.vec,3,3)
		gt.table = cbind(gt.table,allele.a)
		gt.table = cbind(gt.table,allele.b)
	}
	data.list = list("gt.table"=gt.table,"chr.table"=chr.table)
	print(paste(chr.vec[j],"total snps:",as.character(nrow(gt.table))))
	return(data.list)	
}
###############################################################################

#data trimming#
#removes snps that aren't biallelic (new version)
biallele.func = function(data.list) {
	print("biallele.func starting")	
	ref.table=data.list[[2]]
	gt.table=data.list[[1]]
	split.vec = strsplit(ref.table[,4],split="")
	split.length = lapply(split.vec,length)
	rm.sites = which(split.length!=1)
	gt.table=gt.table[-c(rm.sites),]
	ref.table=ref.table[-c(rm.sites),]
	print("biallele.func finished")
	data.list = list("gt.table"=gt.table,"ref.table"=ref.table)
	return(data.list)
} 

#missing data:
missing.data.func = function(data.list) {
	print("starting missing.data.func")	
	gt.table=data.list[[1]]
	ref.table=data.list[[2]]
	rm.vec = 0
	for (i in c(1:nrow(gt.table))) {
		if (length(which(gt.table[i,1:20]=="."))>10) {
			rm.vec=c(rm.vec,i)	
			next	
		}
		if (length(which(gt.table[i,25:44]=="."))>10) {
			rm.vec=c(rm.vec,i)	
			next	
		}
		if (length(which(gt.table[i,45:64]=="."))>10) {
			rm.vec=c(rm.vec,i)	
			next	
		}
		if (length(which(gt.table[i,21:24]=="."))==4) { #removes missing outgroup data
			rm.vec=c(rm.vec,i)	
			next	
		}
		if (length(which(gt.table[i,21:24]=="."))==0) {
			if (length(which(gt.table[i,21:24]=="1"))==2) { #removes if freq is .5
				rm.vec=c(rm.vec,i)	
				next	
			}
		}
		if (length(which(gt.table[i,21:24]=="."))==2) {
			if (length(which(gt.table[i,21:24]=="1"))==1) {	#removes if freq is .5	
				rm.vec=c(rm.vec,i)	
				next	
			}
		}
	}
	ref.table.filt=ref.table[-c(rm.vec),]
	gt.table.filt=gt.table[-c(rm.vec),]
	print(paste(chr.vec[j],"usable snps:",as.character(nrow(gt.table.filt))))
	data.list = list("gt.table.filt"=gt.table.filt,"ref.table.filt"=ref.table.filt)
	return(data.list)
}	

#frequencies
freq.func = function(data.list) {
	#recover()
	print("starting freq.func")	
	gt.table.filt=data.list[[1]]
	ref.table.filt=data.list[[2]]
	og.table=gt.table.filt[,c(21:24)]
	new.og.table = matrix(nrow=nrow(og.table),ncol=ncol(og.table))
	for (z in 1:ncol(og.table)) {
		new.og.table[,z]=suppressWarnings(as.numeric(og.table[,z]))
	}
	og.table=new.og.table
	og.freqs = apply(og.table,1,mean,na.rm=TRUE)
	#allele frequencies for each species
	cyd.freqs = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		sites = as.numeric(gt.table.filt[i,which(gt.table.filt[i,c(1:20)]!=".")])
		cyd.freqs[i] = sum(sites)/length(sites)
	}
	pac.freqs = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		sites = as.numeric(gt.table.filt[i,(which(gt.table.filt[i,c(45:64)]!=".")+44)])
		pac.freqs[i] = sum(sites)/length(sites)
	}
	mel.freqs = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		sites = as.numeric(gt.table.filt[i,(which(gt.table.filt[i,c(25:44)]!=".")+24)])
		mel.freqs[i] = sum(sites)/length(sites)
	}	
	#table building
	og.cons = round(og.freqs)
	og.pol = og.freqs
	og.pol[intersect(which(og.pol!=0),which(og.pol!=1))] = .75
	og.pol[which(og.pol!=.75)] = 1
	#
	freq.table = cbind(ref.table.filt,og.cons)
	freq.table = cbind(freq.table,cyd.freqs)
	freq.table = cbind(freq.table,pac.freqs)
	freq.table = cbind(freq.table,mel.freqs)
	#polarizing 
	freq.table.pol=freq.table
	freq.table.pol[,c(6:8)]= abs(freq.table[,5]-freq.table[,c(6:8)])
	freq.table.pol[,5] = og.pol	
	freq.table[,5] = og.freqs
	colnames(freq.table) = c("chr","pos","ref","alt","og.freqs","cyd.freqs","pac.freqs","mel.freqs")
	colnames(freq.table.pol) = c("chr","pos","ref","alt","og.pol","cyd.freqs","pac.freqs","mel.freqs")
	save(freq.table,file=paste(chr.vec[j], "freqs_raw.RDATA",sep=""))
	save(freq.table.pol,file=paste(chr.vec[j], "freqs_pol.RDATA",sep=""))
	print(paste(chr.vec[j],"finished"))
}

#Get sample sizes
n.func = function(data.list) {
#	recover()
	print("starting n.func")	
	gt.table.filt=data.list[[1]]
	cyd.n = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		cyd.n[i] = length(gt.table.filt[i,which(gt.table.filt[i,c(1:20)]!=".")])
	}
	pac.n = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		pac.n[i] = length(gt.table.filt[i,(which(gt.table.filt[i,c(45:64)]!=".")+44)])
	}
	mel.n = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		mel.n[i] = length(gt.table.filt[i,(which(gt.table.filt[i,c(25:44)]!=".")+24)])
	}	
	n.table = cbind(cyd.n,pac.n)
	n.table = cbind(n.table,mel.n)
	colnames(n.table) = c("cyd.n","pac.n","mel.n")
	save(n.table,file=paste(chr.vec[j], "n_table.RDATA",sep=""))
	
}

#Get heterozygotes
het.func = function(data.list) {
#	recover()
	print("starting het.func")	
	gt.table.filt=data.list[[1]]
	cyd.het = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		gen.vec = gt.table.filt[i,which(gt.table.filt[i,c(1:20)]!=".")]
		length.vec = length(gen.vec)
		evens = seq(2,length.vec,2)
		odds = seq(1,(length.vec-1),2)
		genotypes = cbind(gen.vec[odds],gen.vec[evens])
		hets = length(which((genotypes[,1]==genotypes[,2])==FALSE))
		cyd.het[i] = hets/nrow(genotypes)
	}
	pac.het = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		gen.vec = gt.table.filt[i,(which(gt.table.filt[i,c(45:64)]!=".")+44)]
		length.vec = length(gen.vec)
		evens = seq(2,length.vec,2)
		odds = seq(1,(length.vec-1),2)
		genotypes = cbind(gen.vec[odds],gen.vec[evens])
		hets = length(which((genotypes[,1]==genotypes[,2])==FALSE))
		pac.het[i] = hets/nrow(genotypes)
	}
	mel.het = c(1:nrow(gt.table.filt))
	for (i in c(1:nrow(gt.table.filt))) {
		gen.vec = gt.table.filt[i,(which(gt.table.filt[i,c(25:44)]!=".")+24)]
		length.vec = length(gen.vec)
		evens = seq(2,length.vec,2)
		odds = seq(1,(length.vec-1),2)
		genotypes = cbind(gen.vec[odds],gen.vec[evens])
		hets = length(which((genotypes[,1]==genotypes[,2])==FALSE))
		mel.het[i] = hets/nrow(genotypes)
	}	
	het.table = cbind(cyd.het,pac.het)
	het.table = cbind(het.table,mel.het)
	colnames(het.table) = c("cyd.het","pac.het","mel.het")
	save(het.table,file=paste(chr.vec[j], "het_table.RDATA",sep=""))
}















