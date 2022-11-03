blah <- read.table("/Users/brian/test/1026122_23153_0_0.final.tsv",
                                    sep = "\t", comment.char = "", quote = "", header = T)
x=blah
MUTS=x


which(MUTS$ch_my_pd>=1)
which(is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes))


# bad is n.HGVS is NA
MUTS$ch_my_pd[intersect(which(MUTS$ch_my_pd>=1), which(is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes)))]
ifelse((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes, 1, MUTS$ch_my_pd)[intersect(which(MUTS$ch_my_pd>=1), which(is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes)))]


# okay if n.HGVS not NA
MUTS$ch_my_pd[intersect(which(MUTS$ch_my_pd>=1), which(!is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes)))]
ifelse((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes, 1, MUTS$ch_my_pd)[intersect(which(MUTS$ch_my_pd>=1), which(!is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes)))]


# fix
MUTS$n.HGVSp <- fillna(MUTS$n.HGVSp, 0)
MUTS$n.HGVSc <- fillna(MUTS$n.HGVSc, 0)
MUTS$ch_my_pd[intersect(which(MUTS$ch_my_pd>=1), which(!is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes)))]
ifelse((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes, 1, MUTS$ch_my_pd)[intersect(which(MUTS$ch_my_pd>=1), which(!is.na((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes)))]



