
library(data.table)
library(coloc)
library(reshape)

sbp<-fread('SBP_female_male_bolt_imputed_result_MAF0.01_short.txt',header=T)
dbp<-fread('DBP_female_male_bolt_imputed_result_MAF0.01_short.txt',header=T)
pp<-fread('PP_female_male_bolt_imputed_result_MAF0.01_short.txt',header=T)

top.sbp=read.csv(paste('SBP_topsnp_reigon_annot.csv',sep=''),header=T)
top.dbp=read.csv(paste('DBP_topsnp_reigon_annot.csv',sep=''),header=T)
top.pp=read.csv(paste('PP_topsnp_reigon_annot.csv',sep=''),header=T)

#Example usage
#cvd.data<-fread('example_data_fmd.txt',header=T)
#bp.sex.colocalization(cvd.data,bp.trait='PP',cvd.trait='FMD',size=8656,p=0.3,Type='b',poster.p='H4',gene.pull.method='max',wd=250000, diff=0.5, cutoff=0.5)

bp.sex.colocalization<-function(cvd.data,bp.trait='PP',cvd.trait='CVD',size=1000,p=0.3,Type='b',poster.p='H4',gene.pull.method='mean',wd=250000, diff=0.5, cutoff=0.5){
start_time <- Sys.time()
if(bp.trait=='SBP'){data1=sbp;top=top.sbp}
if(bp.trait=='DBP'){data1=dbp;top=top.dbp}
if(bp.trait=='PP'){data1=pp;top=top.pp}
top=top[order(top$region),]

out.all<-NULL
idx=NULL
for(i in 1:nrow(top)){
    print(i)
	locus.chr=top$CHR[i]
	locus.bp=top$BP[i]
	data.subset<-data1[which(data1$CHR==locus.chr&data1$BP>=locus.bp-wd&data1$BP<=locus.bp+wd),]
	data.subset2<-cvd.data[which(cvd.data$CHR==locus.chr&cvd.data$BP>=locus.bp-wd&cvd.data$BP<=locus.bp+wd),]

	if(nrow(data.subset)>0&nrow(data.subset2)>0){
		bp.male<-data.subset[,c(1:7,13)]
		bp.female<-data.subset[,c(1:2,8:13)]
		colnames(bp.male)<-colnames(bp.female)<-c('chr','BP','rsid','beta','se','pval','maf','SNPID')
		bp.male=bp.male[order(bp.male$pval),]; bp.male=bp.male[!duplicated(bp.male$SNPID),]
		bp.female=bp.female[order(bp.female$pval),];bp.female=bp.female[!duplicated(bp.female$SNPID),]
		input1 <- merge(bp.female, data.subset2, by="SNPID", all=FALSE, suffixes=c("_gwas_female_trait1","_gwas_trait2"))
		input2 <- merge(bp.male, data.subset2, by="SNPID", all=FALSE, suffixes=c("_gwas_male_trait1","_gwas_trait2"))
		m=which(is.na(as.numeric(input1$pval_gwas_trait2))==T)
		if(length(m)>0){input1=input1[-m,];input2=input2[-m,]}
		if(nrow(input1)!=0){
			if(Type=='b'){
				result1 <- coloc.abf(dataset1=list(snp=input1$SNPID, pvalues=input1$pval_gwas_female_trait1,type="quant", N=174664, beta=input1$beta,verbeta=input1$se^2), dataset2=list(snp=input1$SNPID, pvalues=as.numeric(input1$pval_gwas_trait2), beta=input1$BETA,verbeta=input1$SE^2, type="cc", N=size, s=p), MAF=input1$maf)
				result2 <- coloc.abf(dataset1=list(snp=input2$SNPID, pvalues=input2$pval_gwas_male_trait1,  type="quant", N=174664, beta=input2$beta,verbeta=input2$se^2), dataset2=list(snp=input2$SNPID, pvalues=as.numeric(input2$pval_gwas_trait2), beta=input2$BETA,verbeta=input2$SE^2, type="cc", N=size, s=p), MAF=input2$maf)
			}
			if(Type=='q'){
				result1 <- coloc.abf(dataset1=list(snp=input1$SNPID, pvalues=input1$pval_gwas_female_trait1,type="quant", N=174664, beta=input1$beta,verbeta=input1$se^2), dataset2=list(snp=input1$SNPID, pvalues=as.numeric(input1$pval_gwas_trait2), beta=as.numeric(input1$BETA),verbeta=as.numeric(input1$SE)^2, type="quant", N=size), MAF=input1$maf)
				result2 <- coloc.abf(dataset1=list(snp=input2$SNPID, pvalues=input2$pval_gwas_male_trait1,type="quant", N=174664, beta=input2$beta,verbeta=input2$se^2), dataset2=list(snp=input2$SNPID, pvalues=as.numeric(input2$pval_gwas_trait2), beta=as.numeric(input2$BETA),verbeta=as.numeric(input2$SE)^2, type="quant", N=size), MAF=input2$maf)
			}

			out<-c(result1$summary,result2$summary[-1])
			out.all<-rbind(out.all,out)
			idx=c(idx,i)
		}
	}
}
out.all1<-cbind(top[idx,],out.all)
colnames(out.all1)[9:13]=paste('BPF.',colnames(out.all1)[9:13],sep='')
colnames(out.all1)[14:18]=paste('BPM.',colnames(out.all1)[14:18],sep='')
out.all1$BPF.maxPP.H3.H4.abf=apply(cbind(out.all1$BPF.PP.H3.abf,out.all1$BPF.PP.H4.abf),1,max)
out.all1$BPM.maxPP.H3.H4.abf=apply(cbind(out.all1$BPM.PP.H3.abf,out.all1$BPM.PP.H4.abf),1,max)
out.all1$ST=out.all1$BP-wd
out.all1$ED=out.all1$BP+wd

if(poster.p=='H4'){

	if(gene.pull.method=='mean'){
		out.matrix.female<-as.matrix(unlist(by(out.all1$BPF.PP.H4.abf,out.all1$closegene,function(m){mean(m,na.rm=T)})))
		out.matrix.male<-as.matrix(unlist(by(out.all1$BPM.PP.H4.abf, out.all1$closegene,function(m){mean(m,na.rm=T)})))
		out.matrix.gene.st<-as.matrix(unlist(by(out.all1$ST, out.all1$closegene,function(m){min(m,na.rm=T)})))
        out.matrix.gene.ed<-as.matrix(unlist(by(out.all1$ED, out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.gene.chr<-as.matrix(unlist(by(out.all1$CHR, out.all1$closegene,function(m){min(m,na.rm=T)})))
	}

	if(gene.pull.method=='max'){
		out.matrix.female<-as.matrix(unlist(by(out.all1$BPF.PP.H4.abf,out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.male<-as.matrix(unlist(by(out.all1$BPM.PP.H4.abf, out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.gene.st<-as.matrix(unlist(by(out.all1$ST, out.all1$closegene,function(m){min(m,na.rm=T)})))
        out.matrix.gene.ed<-as.matrix(unlist(by(out.all1$ED, out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.gene.chr<-as.matrix(unlist(by(out.all1$CHR, out.all1$closegene,function(m){min(m,na.rm=T)})))
        
	}
}

if(poster.p=='max.H3.H4'){

	if(gene.pull.method=='mean'){
		out.matrix.female<-as.matrix(unlist(by(out.all1$BPF.maxPP.H3.H4.abf,out.all1$closegene,function(m){mean(m,na.rm=T)})))
		out.matrix.male<-as.matrix(unlist(by(out.all1$BPM.maxPP.H3.H4.abf, out.all1$closegene,function(m){mean(m,na.rm=T)})))
		out.matrix.gene.st<-as.matrix(unlist(by(out.all1$ST, out.all1$closegene,function(m){min(m,na.rm=T)})))
        out.matrix.gene.ed<-as.matrix(unlist(by(out.all1$ED, out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.gene.chr<-as.matrix(unlist(by(out.all1$CHR, out.all1$closegene,function(m){min(m,na.rm=T)})))
	}

	if(gene.pull.method=='max'){
		out.matrix.female<-as.matrix(unlist(by(out.all1$BPF.maxPP.H3.H4.abf,out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.male<-as.matrix(unlist(by(out.all1$BPM.maxPP.H3.H4.abf, out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.gene.st<-as.matrix(unlist(by(out.all1$ST, out.all1$closegene,function(m){min(m,na.rm=T)})))
        out.matrix.gene.ed<-as.matrix(unlist(by(out.all1$ED, out.all1$closegene,function(m){max(m,na.rm=T)})))
		out.matrix.gene.chr<-as.matrix(unlist(by(out.all1$CHR, out.all1$closegene,function(m){min(m,na.rm=T)})))
	}
}

print('finalize the result and write to files')

m=which((abs(out.matrix.female[,1]-out.matrix.male[,1])>diff)&(out.matrix.female[,1]>=cutoff&out.matrix.male[,1]<cutoff)|(out.matrix.female[,1]<cutoff&out.matrix.male[,1]>=cutoff))
if(length(m)>0){
	print(paste('Gene regions with suspected sex-specific pleiotropy:', paste(rownames(out.matrix.female)[m],collapse=",")))
	sex.gene.matrix=cbind(out.matrix.gene.chr[m,],out.matrix.gene.st[m,],out.matrix.gene.ed[m,],out.matrix.female[m,],out.matrix.male[m,])
	colnames(sex.gene.matrix)=c('CHR','ST','ED','PP.Female','PP.Male')
	write.csv(sex.gene.matrix,paste(bp.trait,'_',cvd.trait,'_GWAS_colocalization_topSEXGene.csv',sep=''))

}

all.gene.matrix=cbind(out.matrix.gene.chr,out.matrix.gene.st,out.matrix.gene.ed,out.matrix.female,out.matrix.male)
colnames(all.gene.matrix)=c('CHR','ST','ED','PP.Female','PP.Male')
write.csv(all.gene.matrix,paste(bp.trait,'_',cvd.trait,'_GWAS_colocalization.csv',sep=''))

end_time <- Sys.time()
print(end_time - start_time)
return(sex.gene.matrix)
}



plot.region<-function(data,ch,bp,p,locus,pos.chr,pos.st,pos.ed,maxp,alpha){
###locus plot for the entire region
	par(mar=c(5,10,5,2),font.axis=2,font=2,font.lab=2)
	par(mgp=c(1.6,0.5,0))
	subdata<-subset(data,data[,ch]==pos.chr&data[,bp]>=pos.st&data[,bp]<=pos.ed)
	data.subset2<-cvd.data[which(cvd.data$CHR==locus.chr&cvd.data$BP>=locus.bp-wd&cvd.data$BP<=locus.bp+wd),]

	m<-which(subdata[,p]<alpha)
	#maxp<-max(-log10(subdata[,p]),na.rm=T)
	#locus figure	
	plot(subdata[,bp],-log10(subdata[,p]),pch=19,ylab='-log10(P-value)',cex=1.4,cex.lab=1.5,cex.axis=1.2,cex.main=1.8,col='gray31',ylim=c(0,maxp+2),xlim=c(pos.st,pos.ed),xlab='pos(bp.hg19)',main=locus)
	points(subdata[m,bp],-log10(subdata[m,p]),col='brown',pch=19,cex=1.6)
}

##annotation
annot_plot<-function(pos.chr,pos.st,pos.ed){
	size<-pos.ed-pos.st
	par(mar=c(2,10,0,2),font.axis=2,font=2,font.lab=2)
	plot(size,2,xlim=c(0,size),ylim=c(0,7),ann='F',axes='F',type='n',col='white')
	#text(c(-3,-3),c(6.75,5.75),labels=c('Human EST','Human_mRNA'),col='brown',cex=1.1,font=4,srt = 0,adj=1.34,xpd=T)	
	text(c(-3,-3,-3),c(2.75),labels=c('RefGene'),col=c('dark blue'),cex=1.5,font=4,srt = 0,adj=1.85,xpd=T)	

	##DNA annotation
	refgene<-read.table('~/Work/Reference_database/annotation/refGene_hg19.txt',header=F)

	##refgene
	subref<-subset(refgene,(refgene[,3]==paste('chr',pos.chr,sep='')&(refgene[,5]>=pos.st&refgene[,5]<=pos.ed))|(refgene[,3]==paste('chr',pos.chr,sep='')&(refgene[,6]>=pos.st&refgene[,6]<=pos.ed)))
	d<-which(duplicated(subref[,13])=='TRUE')
	if(length(d)!=0){subref[d,13]<-''}
	if(nrow(subref)!=0){
		geneout<-NULL
		#gene figure
		for(i in 1:nrow(subref)){
		tmpstart<-as.character(subref[i,10])
		tmpend<-as.character(subref[i,11])
		start<-unlist(strsplit(tmpstart,',',fixed=T))
		end<-unlist(strsplit(tmpend,',',fixed=T))
		  for(j in 1:length(start)){
			  start1<-as.numeric(start[j])
			  end1<-as.numeric(end[j])
			  rect(start1-pos.st,2.5,end1-pos.st,3.0,col='dark blue',border='dark blue')
			  }	
		lines(c(subref[i,5]-pos.st,subref[i,6]-pos.st),c(2.75,2.75),col='dark blue',lwd=1.6)
		text((subref[i,5]-pos.st+subref[i,6]-pos.st)/2,1,labels=subref[i,13],col='blue',cex=1.5,font=4,srt=90,xpd=T)	  
		if(is.na(subref[i,13])!='TRUE'){geneout<-paste(geneout,',',subref[i,13],sep='')}
			}
		geneout<-sub(',','',geneout)
	}else{geneout<-'NA'}
}


bp.sex.region.locus.plot<-function(outname,cvd.data,bp.trait='PP',cvd.trait='CVD',pos.chr,pos.st,pos.ed){
	if(bp.trait=='SBP'){data1=sbp}
	if(bp.trait=='DBP'){data1=dbp}
	if(bp.trait=='PP'){data1=pp}
	sub1<-subset(data1,CHR==pos.chr&BP>=pos.st&BP<=pos.ed)
	set2<-cvd.data[which(cvd.data$CHR==pos.chr&cvd.data$BP>=pos.st&cvd.data$BP<=pos.ed),]

	locus1=paste(bp.trait,' female-only GWAS',sep='')
	locus2=paste(bp.trait,' male-only GWAS',sep='')
	locus3=paste(cvd.trait,' GWAS',sep='')
	maxp=max(c(-log10(sub1[,6]),-log10(sub1[,11]),-log10(set2[,5]),na.rm=T))
	png(paste(outname,'gene_locus.png',sep=''),width=1000,height=1000)
	layout(matrix(data=c(1,2,3,4), nrow=4, ncol=1), heights=c(3,3,3,3))
	plot.region(sub1,1,2,11,locus1,pos.chr,pos.st,pos.ed,maxp,5e-8)
	plot.region(sub1,1,2,6,locus2,pos.chr,pos.st,pos.ed,maxp,5e-8)
	plot.region(data.subset2,1,2,31,locus3,pos.chr,pos.st,pos.ed,maxp,5e-8)
	annot_plot(pos.chr,pos.st,pos.ed)
	dev.off()
}
