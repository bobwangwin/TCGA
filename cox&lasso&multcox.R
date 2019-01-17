#library(UCSCXenaTools)

## 本项目主要通过使用UCSCXena的原始数据，将多种数据库进行统一整合并进行相关计算，
## 这里通过对胰腺癌的 TCGA GTEx 等数据库，利用均一化处理的数据进行LncRNA表达以及生存分析。

#读取全部数据
library(data.table)
full<-fread("../TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2", header = F,stringsAsFactors = F, sep = '\t')
#cli2 = XenaPrepare("TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz")
#生成TCGA与GTEx的panc数据list
aa<-grep('(TCGA-2J|TCGA-2L|TCGA-2P|TCGA-3A|TCGA-3E|TCGA-5Q|TCGA-F2|TCGA-FB|TCGA-FQ|TCGA-FZ|TCGA-H6|TCGA-H8|TCGA-HV|TCGA-HZ|TCGA-IB|TCGA-KG|TCGA-L1|TCGA-LB|TCGA-M8|TCGA-OE|TCGA-P9|TCGA-PZ|TCGA-Q3|TCGA-RB|TCGA-RL|TCGA-RV|TCGA-S4|TCGA-SD|TCGA-US|TCGA-WF|TCGA-XD|TCGA-XN|TCGA-YB|TCGA-YH|TCGA-YY|TCGA-Z5|TCGA-Z8|TCGA-ZW)',xcol)
ab<-xcol[aa]
tcga_list<-read.table("tcga_list.txt")
gtex_list<-read.table("GTEx_list.txt")
total_list<- rbind(tcga_list,gtex_list)
write.table(total_list,file = "total_list.txt",row.names = F)
#将数据list颠倒并和total list进行match

full[1,1] = "id"
full_t<-t(full)
list <- read.table("total_list.txt",header = T)
m<-full_t[match(list$V1,full_t),]
m<-as.data.frame(m)
m<-na.omit(m)
col_list<-as.character(unlist(full_t[1,]))
colnames(full_t)<-col_list
colnames(m)<-colnames(full_t)
write.table(m, file = "panc_EXP.txt", row.names = F)

# 读取癌组织全表达，并与lncRNA 或miRNA等的基因LIST进行match
lncRNA_list<-read.table("panc_lncRNA_list.txt",header = T, sep = '\t')
pancEXP<-read.table("panc_EXP.txt",header = F,sep = ' ')
pancEXP_t<-t(pancEXP)
pancEXP_t<-as.data.frame(pancEXP_t)
col_list<-as.character(unlist(pancEXP_t[1,]))
colnames(pancEXP_t)<-col_list
pancEXP_tx<-pancEXP_t[-1,]
m<-pmatch(lncRNA_list$ensg,pancEXP_tx$id) # 模糊比对
m<-na.omit(m)
panc_lncRNA_EXP<-pancEXP_tx[m,]
panc_lncRNA_EXPx<-cbind(lncRNA_list$geneid,panc_lncRNA_EXP) #add gene id
write.table(panc_lncRNA_EXPx,"panc_lncRNA_EXP.txt", row.names = F, quote = F, sep = '\t')


# 使用limma进行差异基因筛选,这里需要注意的是，由于
# 下载的数据已经是normalized，因此这里不能再使用Deseq
# 进行数据分析，应使用limma包进行差异分析(可对等考虑：
# 芯片数据的读取结果即已经是normalized数据)
library(limma)
library(pheatmap)
lncRNA_exp<-read.table("panc_lncRNA_EXP_count_unique.txt", header = T)
exprSet <- lncRNA_exp[,2:ncol(x)]
row.names(exprSet)<-lncRNA_exp[,1]
    #设置样本数
condition=factor(c(rep("T",179),rep("N",169)))
design<-model.matrix(~-1+condition)
colnames(design)<-c("Tumor","Normal")
contranst.matrix<-makeContrasts(contrasts="Tumor-Normal",levels=design)
fit<-lmFit(exprSet,design)
fit1<-contrasts.fit(fit,contranst.matrix)
fit2<-eBayes(fit1)
dif<-topTable(fit2,coef="Tumor-Normal",n=nrow(fit2),adjust="BH") 
genesymbol<-rownames(dif)
dif<-cbind(genesymbol,dif)
write.table(dif,file="panc_lncRNA_Foldchange.txt",sep='\t',quote=F,row.names=F)

##### 进行单因素或多因素回归分析

    #将表达矩阵整理成患者行，基因列
expr<-read.table("panc_lncRNA_EXP_count_unique.txt",header = F, sep = '\t')
exprSet<-t(expr)
exprSet<-as.data.frame(exprSet)
xx<-as.character(unlist(exprSet[1,]))
xx[1] = "id"
colnames(exprSet)<-xx
exprSet<-exprSet[-1,]

    #将生存时间列表id与表达矩阵进行匹配
list<- read.table("survival_time.txt",header = T)
list_sur<-exprSet[match(list$id,exprSet$id),]
list_surv<-merge(list,list_sur,by = "id")
write.table(list_surv,"01_panclncRNA_sur_exp.txt", quote = F,row.names = F, sep = '\t') # 手动去除生存期为0天的患者信息

    #进行unicox的计算
library(survival)
inputfile="01_panclncRNA_sur_exp.txt"
lncRNA<-read.table(inputfile,header=T,sep="\t",row.names = 1,check.names = F)

#lncRNAEXP=log2(lncRNA[,3:ncol(lncRNA)]+1)  # 原始数据已经经过log2处理，因此此处不进行
#lncRNA=cbind(lncRNA[,1:2],lncRNAEXP)       # 如果数据为原始count，则要进行log2处理

coxR=data.frame()
coxf<-function(x){
  fmla1 <- as.formula(Surv(survival_time,vital_status)~lncRNA[,x])
  mycox <- coxph(fmla1,data=lncRNA)
}
for(a in colnames(lncRNA[,3:ncol(lncRNA)])){
  mycox=coxf(a)
  coxResult = summary(mycox)
  coxR=rbind(coxR,cbind(lncRNAname=a,HR=coxResult$coefficients[,"exp(coef)"],
                        P=coxResult$coefficients[,"Pr(>|z|)"]))
}
write.table(coxR,"02_coxResult.txt",sep="\t",row.names=F,quote=F)

    # 将p<0.05的基因取出，并进行lasso分析
library(glmnet)
inputfile1="01_panclncRNA_sur_exp.txt" 
inputfile2="02_coxResult005_list.txt" 
inputfile3="survival_time.txt"
data1<-read.table(inputfile1,header = F,sep = "\t",check.names = F)
data2<-read.table(inputfile2,header = T,sep = "\t",check.names = F)
data3<-read.table(inputfile3,header = T,sep = '\t',check.names = F)
data1<-data1[,c(-2,-3)]
data1<-t(data1)
data1<-as.data.frame(data1)
xxx<-as.character(unlist(data1[1,]))
colnames(data1)<-xxx
data1<-data1[-1,]
match_data<-data1[match(data2$id,data1$id),]
write.table(match_data,"03_lassodata.txt",sep = "\t",row.names = F,quote = F) # excel 手动倒置一下
    #将生存时间整合进入lassodata
merger_data<-read.table("03_lassodata.txt",header = T,sep = '\t')
merger_data<-merge(data3,merger_data,by = 'id')
write.table(merger_data,"03_lassoinput.txt",sep = "\t",row.names = F,quote = F)
    #进行lasso分析
inputfile="03_lassoinput.txt"
lncRNA<-read.table(inputfile,header=T,sep="\t",row.names = 1,check.names = F,stringsAsFactors = F) 

#lncRNAEXP=log2(lncRNA[,3:ncol(lncRNA)]+1)
#lncRNA=cbind(lncRNA[,1:2],lncRNAEXP)
lncRNA[,"survival_time"]=lncRNA[,"survival_time"]/365
v1<-as.matrix(lncRNA[,c(3:ncol(lncRNA))])
v2 <- as.matrix(Surv(lncRNA$survival_time,lncRNA$vital_status))
myfit <- glmnet(v1, v2, family = "cox")
pdf("lambda.pdf")
plot(myfit, xvar = "lambda", label = TRUE)
dev.off()

myfit2 <- cv.glmnet(v1, v2, family="cox")
pdf("min.pdf")
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed")
dev.off()

myfit2$lambda.min
coe <- coef(myfit, s = myfit2$lambda.min)
act_index <- which(coe != 0)
act_coe <- coe[act_index]
row.names(coe)[act_index]
write.table(row.names(coe)[act_index],"03lasso_keygene.txt", quote = F,row.names = F)

###  进行多因素回归分析
  #读取 lasso_keygene 信息
inputfile="04mulcoxData.txt" #手动将keygene的数据提取出来并保存 
lncRNA<-read.table(inputfile,header=T,sep="\t",row.names = 1,check.names = F,stringsAsFactors = F) 

#lncRNAEXP=log2(lncRNA[,3:ncol(lncRNA)]+1)
#lncRNA=cbind(lncRNA[,1:2],lncRNAEXP)
lncRNA[,"survival_time"]=lncRNA[,"survival_time"]/365
fmla1 <- as.formula(Surv(survival_time,vital_status)~.)
mycox <- coxph(fmla1,data=lncRNA)
risk_score<-predict(mycox,type="risk",newdata=lncRNA)
risk_level<-as.factor(ifelse(risk_score>median(risk_score),"High","Low"))
write.table(cbind(id=rownames(cbind(lncRNA[,1:2],risk_score,risk_level)),cbind(lncRNA[,1:2],risk_score,risk_level)),"04_risk_score.txt",sep="\t",quote=F,row.names=F)

summary(mycox) # 另外保存

library(survminer)

pdf("forest.pdf")
ggforest(mycox,fontsize = 1)
dev.off()

####  计算模型的 ROC
library(timeROC)
lncRNA<-read.table("04_risk_score.txt",header=T,sep="\t")
predict_3_year<- 3
predict_5_year<- 5

ROC<-timeROC(T=lncRNA$survival_time,delta=lncRNA$vital_status,
             marker=lncRNA$risk_score,cause=1,
             weighting="marginal",
             times=c(predict_3_year,predict_5_year),ROC=TRUE)

pdf("ROC.pdf")
plot(ROC,time=predict_3_year,title=FALSE,lwd=3)
plot(ROC,time=predict_5_year,col="blue",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("AUC of 3 year survival: ",round(ROC$AUC[1],3)),
         paste("AUC of 5 year survival: ",round(ROC$AUC[2],3))),col=c("red","blue"),lwd=3)
dev.off()
