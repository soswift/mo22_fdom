# process each run using R script
setwd("~/Documents/Bioinformatics/Projects/NelsonLabMisc/SGD_fDOM/")
source("R/process_aqualog_functions.R")
source("R/eem_plot_functions.R")


# Process individual runs -------------
# round 1
run1 <- process_aqualog(data_directory ="cabral",
                        run_name = "cabral",
                        sample_key_file = "sgd_run_fdom_cabral.tsv",
                        org_file = "cabral_sample_sheet.csv")

run2 <- process_aqualog(data_directory = "varari",
                        run_name = "varari",
                        sample_key_file = "sgd_run_fdom_varari.tsv")



#  Plot Cabral eems -------------------
c_eem = readRDS("cabral/processed_data/cabral_EEM.rds")

plot_eem(c_eem[[1]], rows_as_names = T)

c_ps = lapply(1:length(c_eem), function(i){
  plot_eem(c_eem[[i]], sample_name = names(c_eem)[[i]], rows_as_names = T)
  
  })

pdf("cabral_eem_plots.pdf")
print(c_ps)
dev.off()


# Plot Cabral blanks
c_blank = read_eems("cabral/processed_data/blanks/")

c_bps = lapply(1:length(c_blank), function(i){
  plot_eem(c_blank[[i]], sample_name = c(1:20)[i], rows_as_names = F)
  
})

pdf("cabral_blank_plots.pdf")
print(c_bps)
dev.off()




# Check Cabral ordination --------------------------
library(vegan)
c_ind = fread("cabral/processed_data/cabral_fDOM_indices_out.csv")



highlight = c(
  "C05_L_D_3_30_0700",
  "C11_L_D_3_30_0700",
  "C16_L_N_3_30_1900",
  "C20_H_N_3_31_0100",
  "C19_L_N_3_30_1900",
  "C20_H_D_3_30_1300",
  "C05_H_N_3_31_0100",
  "C08_L_N_3_30_1900"
)

# organize variables
index_columns = c(
  "CobleA",
  "CobleT",
  "Stedmon_D",
  "Lignin"
)

ordinate_dt = function(dt = c_ind,
                       cols = index_columns){
  set.seed("2022")
  
  ord_dat = as.matrix(dt[ , ..cols])
  row.names(ord_dat) = dt$UniqueID
  
  # ordinate
  raw_ord = metaMDS(ord_dat)
  
  ord_points = as.data.table(raw_ord$points, keep.rownames = T)
  setnames(ord_points, "rn", "UniqueID")
  
  return(ord_points)
  
}

c_ord = ordinate_dt()
c_ord[ ,group_12 := ifelse(UniqueID %in% highlight, "Yes","No")]
c_ord[ , date := gsub(".+(._.._....)","\\1",UniqueID)]


p = ggplot(c_ord, aes(x = MDS1, y = MDS2, label = UniqueID, color = date, shape = group_12))+
  geom_point(size = 5)+
  geom_text(nudge_x = 0.2, alpha = 0.7)

ggsave(plot =p, "p.pdf", width = 10, height = 10, scale = 2)



# check updated data -----------------------

dat = fread("SummaryDataMatrix.csv")

setnames(dat, "Sample Name", "UniqueID")
dat = dat[!grepl("Blank", UniqueID)]

dat = dat[  ,.(UniqueID, `Tyrosine-like`, `Tryptophan-like`, `Marine Humic-like`, `Lignin-like`)]

index_columns = c(
  "Tyrosine-like",
  "Tryptophan-like",
  "Marine Humic-like",
  "Lignin-like"
)

up_ord = ordinate_dt(dat)

up_ord[, date := gsub(".+(._.._....)","\\1",UniqueID)]
up_ord[, site := substr(UniqueID, 1, 1)]

c_p = ggplot(up_ord[site == "C"], aes(x = MDS1, y = MDS2, label = UniqueID, color = date))+
  geom_point(size = 5)+
  geom_text(nudge_x = 0.2, alpha = 0.7)

c_p

v_p = ggplot(up_ord[site == "V"], aes(x = MDS1, y = MDS2, label = UniqueID, color = date))+
  geom_point(size = 5)+
  geom_text(nudge_x = 0.2, alpha = 0.7)
v_p
