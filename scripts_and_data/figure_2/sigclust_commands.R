library("sigclust")
raw = read.table("differential_pathways.tab", header=T, row.names=1, sep="\t")
norm = t(apply(raw, 1, function(x)(x-min(x))/(max(x)-min(x))))

df = t(norm)
df <- df[ order(row.names(df)), ]

# is 2016-02 different?:
labels <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1)
sigclust(df, nsim=1000, labflag=1, label=labels, icovest=3)


# is 2017 the same as 2014/2015?:
df <- df[!grepl("2016",row.names(df)), !grepl("2016", colnames(df))]
labels <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2)

sigclust(df, nsim=1000, labflag=1, label=labels, icovest=3)


