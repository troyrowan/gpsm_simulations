readgwas <- function(filepath){
    #read_tsv(filepath, col_names = TRUE, col_types=cols(rs = col_character())) %>%
    read_tsv(filepath,
             col_names = TRUE) %>%
    mutate(q = qvalue(p_score)$qvalues) %>%

    mutate(pos = ps) %>%
    mutate(rs = paste(chr, pos, sep = ":")) %>%
    select(rs, chr, pos , af, p_score, q)
}


ggqq <- function(pvector){
  pvector = pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
  pdf = data.frame(observed = -log10(sort(pvector,decreasing=FALSE)), expected = -log10(ppoints(length(pvector))))
  qqplotted = ggplot(pdf, aes(expected, observed))+
    geom_point()+
    geom_abline(intercept = 0,
                slope = 1,
                colour = "red")+
    labs(x = expression(paste("Expected ", -log10, "(", italic('p'), ")")),
         y = expression(paste("Expected ", -log10, "(", italic('p'), ")")))
  return(qqplotted)
}

ggmanhattan <- function(inputfile,
         prune = 1.0,
         value = p,
         alpha = 0.5,
         pcol = "p_score",
         colors = c("gray12", "gray55"), 
         sigsnps = NULL,
         sigsnps2 = NULL,
         sigsnpcolor = "red",
         sigsnpcolor2 = "blue"){
  require(qvalue)
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  #Allows p or q value to be interpreted as column name in ggplot call
  v = enexpr(value)
  gwas = inputfile %>% #reads in input file
    select(rs, chr, pos, p = pcol) %>%
    #mutate(chr = case_when(is.null(chr) ~ as.integer(str_split_fixed(snp, ":", n = 2)[,1]),
    #                       TRUE ~ chr)) %>%
    #mutate(pos = case_when(is.null(pos) ~ as.integer(str_split_fixed(snp, ":", n = 2)[,2]),
    #                       TRUE ~ pos)) %>%
    #mutate(chr = as.numeric(str_split_fixed(rs, ":", n = 2)[,1]))%>% #extracts information from SNP name which should be chr:pos
    #mutate(pos = as.numeric(str_split_fixed(rs, ":", n = 2)[,2])) %>% #This holds as long as this was calculated in my imp pipeline
    mutate(q = qvalue(p)$qvalues) %>% #transforms p-values to q-values
    #ifelse(v == "p", filter(p<prune), filter(q<prune))
    filter(., if (v == "p") p < prune else q < prune)

  don = gwas %>%
    group_by(chr) %>%
    summarise(chr_len=max(pos)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas, ., by=c("chr"="chr")) %>%
    arrange(chr, pos) %>%
    mutate( BPcum=pos+tot)

  axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  #v from above
  gwas_plot = ggplot(don, aes(x=BPcum, y=-log10(!!v))) +
    geom_point(aes(color=as.factor(chr)), alpha=alpha, size=1.3) +
    scale_color_manual(values = rep(colors, 29)) +
    scale_x_continuous(label = axisdf$chr[c(TRUE, FALSE)], breaks= axisdf$center[c(TRUE, FALSE)] ) +
    #scale_y_continuous(expand = c(0, 0.01) ) +
    labs(x = "",
         y = case_when(v == "p" ~ expression(paste(-log10, "(", italic('p'), ")")),
                       v == "q" ~ expression(paste(-log10, "(", italic('q'), ")"))))+
    theme_bw() +
    geom_point(data = subset(don, rs %in% sigsnps), color = sigsnpcolor, size = 1.3) +
    geom_point(data = subset(don, rs %in% sigsnps2),color = sigsnpcolor2, size = 1.3) +
    geom_hline(yintercept = case_when(v == "p" ~ -log10(1e-5),
                                      v == "q" ~ -log10(0.1)),
               color = "red",
               size = 0.5) +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
    )+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )+
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 2.5, face = "bold"))
  return(gwas_plot)
}
