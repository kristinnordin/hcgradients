md <- function(x){
  dens <- density(x)
  peak <- dens$x[which.max(dens$y)]
  return(peak)
}

all_ky<-data.frame()

ky %>% dplyr::select(Visual_L) %>% filter(complete.cases(.)) %>% nrow(.)

for(sheet in c("LG1","RG1","LG2","RG2","LG3","RG3")){
# for(sheet in c("O1_LHC_G1","O1_RHC_G2","O2_LHC_G1","O2_RHC_G2")){
    
  ky<-read.xlsx("Gradients/Hippocampus/Kristin_Paper/Yeo_values221102_young.xlsx",sheet=sheet)

  df<-ky %>% pivot_longer(seq(1,ncol(ky)),names_to="region",values_to="eigen") %>%
    mutate(Hemi=ifelse(grepl("_L",region),"lh","rh")) %>%
    mutate(region=gsub("_L|_R","",region)) %>% #select(region) %>% unique()
    #mutate(region=factor(region,levels=c("Default","Limbic","Somatomotor","DorsalAttention","Visual","Frontoparietal","VentralAttention"))) %>%
    filter(!is.na(eigen)) %>%
    mutate(grad=sheet) #%>%
  

  max_density <- df %>%
    ddply(.(region,grad),summarise,max_density = md(eigen))

  all_ky<-rbind(all_ky,left_join(df,max_density,by=c("grad","region")))
  
}

parula_colors <-read.csv2("Gradients/Hippocampus/Kristin_Paper/hex_and_rgb_v1.1.1/parula.csv")
n<-nrow(parula_colors)
parula_colors[n,1]

# parula_colors[n,1]

unique(all_ky$region)
all_ky<-all_ky %>% mutate(region=gsub("DorsalAttention","DAN",region)) %>%
  mutate(region=gsub("VentralAttention","VAN",region)) %>%
  mutate(region=gsub("Frontoparietal","FPN",region)) %>%
  mutate(region=gsub("Default","DMN",region)) 

pp<-c()
for(g in unique(all_ky$grad)){
  # g<-"O1_LHC_G1"
  p<-all_ky %>%
    filter(grad==g) %>% 
    mutate(HemiG=ifelse(grepl("L",grad),"Left","Right")) %>%
    mutate(grad=gsub("L|R","",grad)) %>%
    # mutate(grad=gsub("_LHC_|_RHC_","",grad)) %>%
    # mutate(gr=factor(ifelse(grepl("O1",grad),"Youth-like","Aged"),levels = c("Youth-like","Aged"))) %>%
    # mutate(grad=gsub("O1|O2","",grad)) %>%
    # mutate(grad=factor(grad)) %>% 
    # group_by(grad) %>%
    arrange(desc(max_density)) %>% #select(region,max_density) %>% unique()
    mutate(region=factor(region,levels=unique(region))) %>%
    # select(region,max_density) %>% 
    ggplot(aes(eigen,region,fill=stat(x)))+
    # geom_density_ridges(aes(fill=ref),scale=15,alpha=0.9)+
    geom_density_ridges_gradient(scale=5)+
    scale_fill_gradient2(low = parula_colors[1,1], mid = parula_colors[floor(n/2),1], high = "#FFFF00" , midpoint = 0.5, limit = c(-0.1, 1.1), space = "Lab", guide = "colourbar")+
    # scale_fill_gradient2(low = "blue", mid = "green", high = "yellow", midpoint = 0.5, limit = c(-0.1, 1.1), space = "Lab", guide = "colourbar")+
    # geom_density_ridges_gradient(scale=5)+scale_fill_viridis_c(option="C")+
    facet_wrap(~grad+HemiG,ncol = 2,scales = "free")+
    theme_pubr()+theme(axis.line = element_blank(),axis.ticks = element_blank(),legend.position = "none",strip.background = element_blank())+ylab("")+xlab("")

  ggsave(paste0("Gradients/Hippocampus/Kristin_Paper/test/YOUNG_",g,".pdf"),width=10,height = 10,units = "cm")
}

