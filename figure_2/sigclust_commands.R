# for bins:
library("sigclust")
raw = read.table("bin_abundances.tab")

std = 1000*sweep(raw,2,colSums(raw),`/`)
norm = t(apply(std, 1, function(x)(x-min(x))/(max(x)-min(x))))
has_na <- apply(norm, 1, function(x){any(is.na(x))})
norm <- norm[!has_na,]

df = t(norm)
df <- df[ order(row.names(df)), ]
labels <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
sigclust(df, nsim=100, labflag=1, label=labels, icovest=3)



# for contigs:
library("sigclust")
raw <- read.table("long_contigs_5kb.tab")
std <- 1000*sweep(raw,2,colSums(raw),`/`)
norm = t(apply(std, 1, function(x)(x-min(x))/(max(x)-min(x))))
has_na <- apply(norm, 1, function(x){any(is.na(x))})
norm <- norm[!has_na,]

rand <- norm[sample(nrow(norm), 5000), ]
df <- t(rand)
df <- df[ order(row.names(df)), ]

labels <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
sigclust(df, nsim=100, labflag=1, label=labels, icovest=3)


