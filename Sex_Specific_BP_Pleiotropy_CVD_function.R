
library(data.table)
library(coloc)
library(reshape)

sbp.f<-fread('UKB_SBP_female_bolt_imputed_result_MAF0.01.tsv',header=T)
sbp.m<-fread('UKB_SBP_male_bolt_imputed_result_MAF0.01.tsv',header=T)
dbp.f<-fread('UKB_DBP_female_bolt_imputed_result_MAF0.01.tsv',header=T)
dbp.m<-fread('UKB_DBP_male_bolt_imputed_result_MAF0.01.tsv',header=T)
pp.f<-fread('UKB_PP_female_bolt_imputed_result_MAF0.01.tsv',header=T)
pp.m<-fread('UKB_PP_male_bolt_imputed_result_MAF0.01.tsv',header=T)


top.sbp=read.csv(paste('SBP_topsnp_reigon_annot.csv',sep=''),header=T)
top.dbp=read.csv(paste('DBP_topsnp_reigon_annot.csv',sep=''),header=T)
top.pp=read.csv(paste('PP_topsnp_reigon_annot.csv',sep=''),header=T)

bp.sex.colocalization<-function(cvd.data,bp.trait='PP',cvd.trait='CVD',size=1000,p=0.3,Type='b',poster.p='H4',gene.pull.method='mean',wd=250000, diff=0.5, cutoff=0.5){
start_time <- Sys.time()
if(bp.trait=='SBP'){data1.f=sbp.f;data1.m=sbp.m;top=top.sbp}
if(bp.trait=='DBP'){data1.f=dbp.f;data1.m=dbp.m;top=top.dbp}
if(bp.trait=='PP'){data1.f=pp.f;data1.m=pp.m;top=top.pp}
top=top[order(top$region),]

out.all<-NULL
idx=NULL
for(i in 1:nrow(top)){
    print(i)
	locus.chr=top$CHR[i]
	locus.bp=top$BP[i]
	bp.male<-data1.m[which(data1.m$chromosome==locus.chr&data1.m$base_pair_location>=locus.bp-wd&data1.m$base_pair_location<=locus.bp+wd),]
	bp.female<-data1.f[which(data1.f$chromosome==locus.chr&data1.f$base_pair_location>=locus.bp-wd&data1.f$base_pair_location<=locus.bp+wd),]
	bp.male$MAF=bp.male$effect_allele_frequency;bp.male$MAF[bp.male$effect_allele_frequency>0.5]=1-bp.male$effect_allele_frequency[bp.male$effect_allele_frequency>0.5]
	bp.female$MAF=bp.female$effect_allele_frequency;bp.female$MAF[bp.female$effect_allele_frequency>0.5]=1-bp.female$effect_allele_frequency[bp.female$effect_allele_frequency>0.5]
	data.subset2<-cvd.data[which(cvd.data$CHR==locus.chr&cvd.data$BP>=locus.bp-wd&cvd.data$BP<=locus.bp+wd),]

	if(nrow(bp.male)>0&nrow(data.subset2)>0){
		bp.male=bp.male[order(bp.male$p_value),]; bp.male=bp.male[!duplicated(bp.male$pos.hg19),]
		bp.female=bp.female[order(bp.female$p_value),];bp.female=bp.female[!duplicated(bp.female$pos.hg19),]
		input1 <- merge(bp.female, data.subset2, by="pos.hg19", all=FALSE, suffixes=c("_gwas_female_trait1","_gwas_trait2"))
		input2 <- merge(bp.male, data.subset2, by="pos.hg19", all=FALSE, suffixes=c("_gwas_male_trait1","_gwas_trait2"))
		m=which(is.na(as.numeric(input1$pval))==T)
		if(length(m)>0){input1=input1[-m,];input2=input2[-m,]}
		if(nrow(input1)!=0){
			if(Type=='b'){
				result1 <- coloc.abf(dataset1=list(snp=input1$pos.hg19, pvalues=input1$p_value,type="quant", N=174664, beta=input1$beta,verbeta=input1$standard_error^2), dataset2=list(snp=input1$pos.hg19, pvalues=as.numeric(input1$pval), beta=input1$BETA,verbeta=input1$SE^2, type="cc", N=size, s=p), MAF=input1$MAF)
				result2 <- coloc.abf(dataset1=list(snp=input2$pos.hg19, pvalues=input2$p_value,  type="quant", N=174664, beta=input2$beta,verbeta=input2$standard_error^2), dataset2=list(snp=input2$pos.hg19, pvalues=as.numeric(input2$pval), beta=input2$BETA,verbeta=input2$SE^2, type="cc", N=size, s=p), MAF=input2$MAF)
			}
			if(Type=='q'){
				result1 <- coloc.abf(dataset1=list(snp=input1$pos.hg19, pvalues=input1$p_value,type="quant", N=174664, beta=input1$beta,verbeta=input1$standard_error^2), dataset2=list(snp=input1$pos.hg19, pvalues=as.numeric(input1$pval), beta=as.numeric(input1$BETA),verbeta=as.numeric(input1$SE)^2, type="quant", N=size), MAF=input1$MAF)
				result2 <- coloc.abf(dataset1=list(snp=input2$pos.hg19, pvalues=input2$p_value,type="quant", N=174664, beta=input2$beta,verbeta=input2$standard_error^2), dataset2=list(snp=input2$pos.hg19, pvalues=as.numeric(input2$pval), beta=as.numeric(input2$BETA),verbeta=as.numeric(input2$SE)^2, type="quant", N=size), MAF=input2$MAF)
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
	subdata<-data[which(CHR==pos.chr&BP>=pos.st&BP<=pos.ed),]
	#maxp<-max(-log10(subdata[,p]),na.rm=T)
	#locus figure	
	plot(subdata[,bp],-log10(subdata[,p]),pch=19,ylab='-log10(P-value)',cex=1.4,cex.lab=1.5,cex.axis=1.2,cex.main=1.8,col='gray31',ylim=c(0,maxp+2),xlim=c(pos.st,pos.ed),xlab='pos(bp.hg19)',main=locus)
	abline(h=40, col="blue")
}

##annotation
annot_plot<-function(pos.chr,pos.st,pos.ed){
	size<-pos.ed-pos.st
	par(mar=c(2,10,0,2),font.axis=2,font=2,font.lab=2)
	plot(size,2,xlim=c(0,size),ylim=c(0,7),ann='F',axes='F',type='n',col='white')
	#text(c(-3,-3),c(6.75,5.75),labels=c('Human EST','Human_mRNA'),col='brown',cex=1.1,font=4,srt = 0,adj=1.34,xpd=T)	
	text(c(-3,-3,-3),c(2.75),labels=c('RefGene'),col=c('dark blue'),cex=1.5,font=4,srt = 0,adj=1.85,xpd=T)	

	##DNA annotation
	refgene<-read.table('refGene_hg19.txt',header=F)

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


bp.sex.region.locus.plot<-function(outname,cvd.data,bp.trait='PP',cvd.trait='CVD',locus.chr,locus.st,locus.ed){
	print(outname)
	if(bp.trait=='SBP'){data1.f=sbp.f;data1.m=sbp.m;top=top.sbp}
	if(bp.trait=='DBP'){data1.f=dbp.f;data1.m=dbp.m;top=top.dbp}
	if(bp.trait=='PP'){data1.f=pp.f;data1.m=pp.m;top=top.pp}
	bp.male<-data1.m[which(data1.m$chromosome==locus.chr&data1.m$base_pair_location>=locus.st&data1.m$base_pair_location<=locus.ed),]
	bp.female<-data1.f[which(data1.f$chromosome==locus.chr&data1.f$base_pair_location>=locus.st&data1.f$base_pair_location<=locus.ed),]

	sub2<-cvd.data[which(cvd.data$CHR==locus.chr&cvd.data$BP>=locus.st&cvd.data$BP<=locus.ed),]
	locus1=paste(bp.trait,' female-only GWAS',sep='')
	locus2=paste(bp.trait,' male-only GWAS',sep='')
	locus3=paste(cvd.trait,' GWAS',sep='')
	maxp=max(c(-log10(bp.male$p_value),-log10(bp.female$p_value),-log10(sub2$pval)),na.rm=T)
	png(paste(outname,'gene_locus.png',sep=''),width=1000,height=1000)
	layout(matrix(data=c(1,2,3,4), nrow=4, ncol=1), heights=c(3,3,3,3))
	plot(bp.female$base_pair_location,-log10(bp.female$p_value),pch=19,ylab='-log10(P-value)',cex=1.4,cex.lab=1.5,cex.axis=1.2,cex.main=1.8,col='purple',ylim=c(0,maxp+2),xlim=c(locus.st,locus.ed),xlab='pos(bp.hg19)',main=locus1)
	plot(bp.male$base_pair_location,-log10(bp.male$p_value),pch=19,ylab='-log10(P-value)',cex=1.4,cex.lab=1.5,cex.axis=1.2,cex.main=1.8,col='blue',ylim=c(0,maxp+2),xlim=c(locus.st,locus.ed),xlab='pos(bp.hg19)',main=locus2)
	plot(sub2$BP,-log10(sub2$pval),pch=19,ylab='-log10(P-value)',cex=1.4,cex.lab=1.5,cex.axis=1.2,cex.main=1.8,col='black',ylim=c(0,max(-log10(sub2$pval))+1),xlim=c(locus.st,locus.ed),xlab='pos(bp.hg19)',main=locus3)
	
	annot_plot(locus.chr,locus.st,locus.ed)
	dev.off()
}
