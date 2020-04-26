library(XML)
library(RCurl)
library(rlist)
library(hrbrthemes)
library(ggplot2)
library(viridis)
library(plotly)
library(htmlwidgets)

#TODO https://www.ebi.ac.uk/gwas/genes/FOXO3



#------trait IDs to traid description---------
theurl <- getURL("http://geneatlas.roslin.ed.ac.uk/trait/"
                 ,.opts = list(ssl.verifypeer = FALSE) )
trait.http <- strsplit(theurl, split = '\n')[[1]]
trait.http.option <- trait.http[grep('option value', trait.http)]

library(dplyr)
trait_id_to_desc <- 
  trait.http.option %>% 
  gsub(pattern = '    <option value=\"', replacement = '', x = .) %>% 
  gsub(pattern = '</option>', replacement = '', x = .) %>% 
  strsplit(., split = '\">') %>% 
  do.call(rbind.data.frame, .)
trait_id_to_desc[] <- lapply(trait_id_to_desc, as.character)
colnames(trait_id_to_desc) <- c('web-id', 'trait-code')

#-----trait table-----------

traits.table.url <- getURL("http://geneatlas.roslin.ed.ac.uk/traits-table/",
                 .opts = list(ssl.verifypeer = FALSE) )
trait.table <- 
  traits.table.url %>% 
  readHTMLTable(.)

#convert to df
trait.table <- data.frame(trait.table[[1]])

# convert to character  
trait.table[] <- lapply(trait.table, as.character)


trait.table.merged <- merge(trait_id_to_desc, trait.table, by.x = 'trait-code', 
                     by.y = 'Description', all = TRUE)


#--------------Get info for regions and so on------
#TODO Rename it
# Select the region for GRCh37 assembly in kilobases
# foxo start and stop from https://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXO3
# chr_start <- 108881
# chr_stop <- 109006
chr_start <- 108855
chr_stop <- 109037
chr <- 6
traits <- trait.table.merged$`web-id`

# counts to 5 t prevent overloading webpage
N <- 1
# for (i in 1:length(traits)) {
gwas.tables <- list()
j <- 1

for (i in 1:length(traits)) {
  if (N == 6) {
    url <- paste0("http://geneatlas.roslin.ed.ac.uk/region/?traits=", 
                  trait_field,
                  "&esymbol=&gextrawind=50&maxregion=",
                  chr_stop, "&minregion=",
                  chr_start, "&chrom=", chr,
                  "&representation=table")
    theurl <- getURL(url,.opts = list(ssl.verifypeer = FALSE) )
    gwas.tables[[j]] <- readHTMLTable(theurl)
    j <- j+1
    Sys.sleep(2)
    N <- 1
  }
  if (N == 1)
    trait_field <- traits[i]
  else 
    trait_field <- paste0(trait_field, "%2C", traits[i])
  N <- N+1
}

# make dataframe from list
gwas.tables.unlisted <- gwas.tables %>% 
  lapply(., function(df) {do.call(rbind.data.frame,df)}) %>% 
  do.call(rbind.data.frame, .)

# create -log10 p-value
gwas.tables.unlisted$log10 <- round(-log10(
  as.numeric(as.character(gwas.tables.unlisted$`p-value`))), digits = 1)
gwas.tables.for_plot <- gwas.tables.unlisted[,c(1,2,4,11)]
gwas.tables.for_plot[,1:3] <- lapply(gwas.tables.for_plot[,1:3], as.character)
gwas.tables.for_plot$Position <- as.numeric(gwas.tables.for_plot$Position)

#create filter on -log10(p.val)
filters <- gwas.tables.for_plot %>% 
  group_by(Trait) %>% 
  summarise(Q=max(log10)) %>% 
  arrange(desc(Q)) %>% 
  filter(Q > 4)
gwas.tables.for_plot <- gwas.tables.for_plot[gwas.tables.for_plot$Trait %in%
                                               filters$Trait,]

filters_imp <- gwas.tables.for_plot %>% 
  group_by(Trait) %>% 
  summarise(Mi=median(`imp. score`)) %>% 
  filter(Mi > 0.95)
gwas.tables.for_plot <- gwas.tables.for_plot[gwas.tables.for_plot$Trait %in% 
                                               filters_imp$Trait, ]

# visualise
library(hrbrthemes)
library(ggplot2)
library(viridis)
p <- ggplot(gwas.tables.for_plot, aes(Variant, Trait, fill=log10)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45))
p

fig <- ggplotly(p)

fig

ggplotly(p, tooltip="text")

# save the widget

saveWidget(p, file="ggplotlyHeatmap.html")

#----only CVD--------------
gwas.tables.for_plot <- gwas.tables.unlisted[,c(1,2,4,11)]
gwas.tables.for_plot[,1:3] <- lapply(gwas.tables.for_plot[,1:3], as.character)
gwas.tables.for_plot$Position <- as.numeric(gwas.tables.for_plot$Position)

#create filter on -log10(p.val)
gwas.tables.for_plot <- gwas.tables.for_plot[grep(pattern = 'heart', gwas.tables.for_plot$Trait), ]

# visualise

p <- ggplot(gwas.tables.for_plot, aes(Variant, Trait, fill=log10)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45))
p

fig <- ggplotly(p)


ggplotly(p, tooltip="text")


#--------parse By significance table-----------
sign.tables <- list()
# 284:344
traits_heart <- trait.table.merged$`web-id`[grep('heart', 
                                                 trait.table.merged$`trait-code`) ]
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I2',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I1',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I3',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I4',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I5',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I6',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I7',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I8',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- c(traits_heart, 
                  trait.table.merged$`web-id`[grep('I9',
                                                   trait.table.merged$`trait-code`) ])
traits_heart <- unique(c(traits_heart, 27, 279))

for (k in 1:length(traits_heart)) {
  sign.url <- paste0('http://geneatlas.roslin.ed.ac.uk/bysignificance/?traits=', 
                     traits_heart[k],
                     '&pvalue=1e-8')
  theurl <- getURL(sign.url,.opts = list(ssl.verifypeer = FALSE) )
  sign.tables[[k]] <- readHTMLTable(theurl)
  Sys.sleep(1)
}

# make dataframe from list
sign.tables.unlisted <- sign.tables %>% 
  lapply(., function(df) {do.call(rbind.data.frame,df)}) %>%
  do.call(rbind.data.frame, .)
sign.tables.unlisted[] <- lapply(sign.tables.unlisted, as.character)
sign.tables.unlisted$pv <- as.numeric(sign.tables.unlisted$pv)

# Nick
library("CMplot")      #install.packages("CMplot")
gwas.tables.for_plot$chr <- "6"
gwas.tables.for_plot$p.value <- 10^(-gwas.tables.for_plot$log10)
gwas.tables.manh <- gwas.tables.for_plot[,c(1,5,2,6)]

sign.tables.for_plot <- sign.tables.unlisted[grep('heart', 
                                                  x = sign.tables.unlisted$Trait), ]

sign.tables.manh <- sign.tables.for_plot[,c(2,3,4,7)]

colnames(gwas.tables.manh) <- colnames(sign.tables.manh)
tables.manh <- rbind(sign.tables.manh, gwas.tables.manh)

# gwas.tables.manh$Variant <- as.factor(gwas.tables.manh$Variant)
# gwas.tables.manh$Trait <- as.factor(gwas.tables.manh$Trait)

SNPs <- unique(gwas.tables.for_plot$Variant)
p1 = CMplot(tables.manh, plot.type="c", multracks=FALSE, type="p", band=2,
            r=1.7, cir.legend=TRUE, highlight.type = 'h',
            col=c("grey30","grey60"),
            # threshold=c(1e-9,1e-5), signal.col=c("red","green"),
            bin.size=1e2, LOG10=TRUE, highlight=SNPs,
            amplify = FLASE, highlight.cex = 1,
            outward=TRUE, cir.legend.col="black", cir.chr.h=1.3 ,
            chr.den.col="black", file="jpg",highlight.col = "green",
            memo="", dpi=600)

p1 = CMplot(gwas.tables.manh, plot.type="c",
            r=1.1, cir.legend=TRUE,bin.size=1e2, LOG10=FALSE, highlight=SNPs,
            outward=TRUE, cir.legend.col="black", cir.chr.h=1.3 ,chr.den.col="black", file="jpg",
            memo="", dpi=300)
p1