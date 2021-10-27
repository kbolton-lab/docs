install.packages("dplyr")
read.table

?binom.test()
?prop.test

?t.test
test <- read.table("~/Downloads/CHAJ-113518-243.indel", header =T)
test2 <- read.table("~/Downloads/CHAJ-1150-304.indel", header =T)
test.snp <- read.table("~/Downloads/CHAJ-113518-243.snp", header =T)
test.snp$var[1]
nchar(test.snp$var)>1
test.snp[nchar(test.snp$var)>1,]


HIV <- read.csv("~/Bolton/HIV/id.csv", header=T)
paste0(HIV$instrument_data_id, collapse=",")
length(unique(HIV$instrument_data_id))

library(pwr)
pwr.p.test()




power.prop.test()
install.packages("pwr")
pwr.2p.test(n=161, sig.level=.05, power=0.8, alternative = "two.sided")
pwr.f2.test()

binom.test(x=0, n=161, p = .05, alternative = "less")
1-0.1982743-0.3224460-0.2605624

dbinom(4,161,.01)
dbinom(6,267,.01)

1-pbinom(5, 161, .01)
1-pbinom(3 , 201000, .01)

qbinom(.95, 161, .01)
# UKBB
qbinom(0.95, 1992, .01)
# UKBB
qbinom(0.95, 4785, .01)
for file in $(ls *.final.tsv); do
  sample=$(basename $file .final.tsv)
  awk -F'\t' '$202 < 1.260958e-09 {print}' OFS='\t' $file > $sample.pon.pass.tsv
done


#total_people <- c(2547, 3423, 4546, 8487, 10355, 12693, 11933, 10534, 8882, 5991, 4136, 1935)
#variant <- c(9, 0, 4, 14, 20, 17, 15, 5, 4, 4, 2, 1)
ages <-  c(15, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 90)

all_ages <- c(rep(32.5, 1428),
            rep(37.5, 1606),
            rep(42.5, 1905),
            rep(47.5, 3228),
            rep(52.5, 4577),
            rep(57.5, 3924),
            rep(62.5, 3656),
            rep(67.5, 3194),
            rep(72.5, 2153),
            rep(77.5, 1283))

[6, 1, 2, 7, 11, 6, 5, 4, 1, 1]

variant_ages <- c(
              rep(32.5, 6),
              rep(37.5, 1),
              rep(42.5, 2),
              rep(47.5, 7),
              rep(52.5, 11),
              rep(57.5, 6),
              rep(62.5, 5),
              rep(67.5, 4),
              rep(72.5, 1),
              rep(77.5, 1)
             )


          
h <- hist(all_ages, breaks = 20, freq = F)
xfit<-seq(min(ages),max(ages),length=40)
yfit<-dnorm(xfit,mean=mean(all_ages),sd=sd(all_ages))
#yfit <- (yfit*diff(h$mids[1:2])*length(all_ages))
lines(xfit, yfit, col="blue", lwd=2)
h$density = h$counts/sum(h$counts)*100
#lines(h,freq=FALSE)

xfit2<-seq(min(ages),max(ages),length=40)
yfit2<-dnorm(xfit2,mean=mean(variant_ages),sd=sd(variant_ages))
lines(xfit2, yfit2, col="red", lwd=2)


## statistically significant
pnorm((mean(variant_ages)-mean(all_ages))/(sd(all_ages)/sqrt(length(variant_ages))))


# H1 = mean(variant_ages) - mean(all_ages) > 0
## H1 says variant ages are older that all ages
t.test(variant_ages, all_ages, alternative = "greater", var.equal = FALSE)
ld discriminant analysis
wilcox.test(variant_ages, all_ages, alternative = "less")



mysql> grant all privileges on test.* to 'brian'@'localhost' identified by 's3kr1t';
