#' Title
#'
#' @param nom 
#' @param x 
#' @param pos 
#' @param Est_cnv 
#' @param Sel 
#'
#' @return
#' @export
#'
#' @examples
graphes <- function(nom, x, pos, Est_cnv, Sel, Final ){
  
data_rep <-   tibble(chr = x, x = pos) %>% 
    ggplot(aes(x = x, y = chr)) + geom_point() + 
    labs(y = nom)
  
densities <- tibble(f0x_est = Est_cnv$Em$f0x, f1x_est = Est_cnv$Em$f1x, chr= x) %>% 
  gather(-chr, key ="Densities", value = "value") %>% 
  ggplot(aes(x = chr, y = value, color =Densities)) + geom_line() 

Tb <- enframe(x) %>% 
  rowid_to_column() %>% 
  mutate(pos = pos)
seuil <- 0.1
var1 <- "'{'"

var2 <- "'}'"
selection <- Sel %>% mutate(chr= list(Tb), 
               chr = map2(Sel, chr, ~ mutate(.y, sel = case_when(rowid %in% .x~"Selected",
                                                                   TRUE~ "V")))) %>% 
  unnest(chr) %>% 
  filter(Nom!="egfr") %>% 
  mutate( Nom = case_when( Nom == "pval_tresh" ~  glue("S(X) == bgroup({var1},p <", seuil,",{var2})"), 
                           Nom == "sel_viter_est" ~ "S(X) == bgroup('{',viterbi == 1,'}')",
                           Nom == "sel_viter_or" ~ "S(X) == bgroup('{',viterbi[or] == 1,'}')",
                           Nom == "sel_viter_min_size" ~ "S(X) == bgroup('{','viterbi = 1 & size > 500','}')",
                           
                           Nom == "H0_and_sel_viter_est" ~ "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')",
                           Nom == "H1" ~  "S(X) == H1", 
                           Nom == "pval_ord" ~  "S(X) == bgroup('{',p < th,'}')",
                           Nom == "lfdr_tresh" ~  glue("S(X) == SC(",seuil,")"),
                           Nom == "lfdr_tresh_pval" ~  "S(X) == SC(FDR[p])",
                           Nom == "lfdr_tresh_viter"~  "S(X) == SC(FDR[v])",
                           Nom == "all" ~ "All"), 
          Nom = factor(Nom, levels = c(glue("S(X) == bgroup({var1},p <", seuil,",{var2})"),
                                       "S(X) == SC(FDR[p])",
                                       glue("S(X) == SC(",seuil,")"), 
                                       "S(X) == bgroup('{',p < th,'}')",
                                       "S(X) == bgroup('{',viterbi == 1,'}')", 
                                       "S(X) == bgroup('{',viterbi[or] == 1,'}')",
                                       "S(X) == SC(FDR[v])",
                                       "S(X) == bgroup('{','viterbi = 1 & size > 500','}')", 
                                       "All"))) %>% 
  ggplot(aes(x = pos, y = value, color = sel)) +
  geom_point(size = 0.2) + 
  facet_grid(Nom~., labeller  = label_parsed) +
  scale_color_manual(values = c("orange", "grey20"), 
                     labels = c("Selected"," Not Selected")) + 
  labs( x= "Position", y = "Statistics", color ="")


 Res <- Final %>% mutate(al = 0.1, Size = Size.x) %>% 
  res_FDR(seuil = 0.1, min_size = 500)
 IC <-Res %>% 
  filter(!is.na(Nom)) %>% 
  mutate(
    V_HMM_boot_qdemi = case_when(
      is.na(V_HMM_boot_qdemi)~ V_HMM_boot_samesel, 
      TRUE ~ V_HMM_boot_qdemi
    ),
    V_HMM_small_boot_qdemi = case_when(
      is.na(V_HMM_small_boot_qdemi)~ V_HMM_small_boot_samesel, 
      TRUE ~ V_HMM_small_boot_qdemi
    )
  ) %>% 
  rename(V_small_simes = V_simes_small) %>% 
  select(Nom, 
         # starts_with("V_simes"), V_small_simes,
         V_HMM_small_boot_qdemi, V_HMM_boot_qdemi,
         V_HMM_small_est, V_HMM_est, FDR_est) %>% 
  gather(-Nom, -FDR_est, -V_HMM_boot_qdemi, 
         -V_HMM_est, 
         # -V_simes,
         value ="high", key ="Type1") %>% 
  gather(-Nom, -FDR_est, -high, -Type1,
         value ="small", key ="Type2") %>% 
  mutate(Type1 = gsub("small_","", Type1)) %>% 
  # mutate(Type2 = gsub("V_", "", Type2)) %>% 
  filter(Type1 == Type2) %>% 
  ggplot(aes(x = Nom, ymin = small,
             ymax = high, color = Type1)) + 
  geom_errorbar() + coord_flip()+
  geom_point(aes(y = FDR_est))+
  scale_x_discrete(labels = function(l)parse(text=l)) +
  labs(x ="", y = "FDP") + 
  scale_color_nejm(
    labels = c("Bootstrap", "Plugin"))
return(list(IC, data_rep, densities, selection))
}