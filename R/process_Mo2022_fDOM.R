# process each run using R script
setwd("~/Documents/Bioinformatics/Projects/NelsonLabMisc/Moorea22_fDOM/")
source("R/src/process_aqualog_functions.R")
source("R/src/eem_plot_functions.R")
library(data.table)



run_dirs = c(
  
  SGD_cabral = "data/projects/SGD_fDOM/cabral/",
  SGD_varari = "data/projects/SGD_fDOM/varari/",
  RRC_1 = "data/projects/RRCoral_fDOM/Run1/",
  RRC_2 = "data/projects/RRCoral_fDOM/Run2/",
  Spiffy = "data/projects/Spiffy_fDOM/spiffy/",
  SL = "data/projects/Spiffy_Lang_fDOM/spiffy_lang/",
  Diel_1= "data/projects/Diel/run1/diel/",
  Diel_2 = "data/projects/Diel/run2/",
  ATI_1 = "data/projects/Around_The_Island_fDOM/around_the_island/run1/",
  ATI_2 =  "data/projects/Around_The_Island_fDOM/around_the_island/run2/",
  ATI_3 =  "data/projects/Around_The_Island_fDOM/around_the_island/run3/",
  ATI_4 =  "data/projects/Around_The_Island_fDOM/around_the_island/run4/",
  ATI_5 =  "data/projects/Around_The_Island_fDOM/around_the_island/run5/",
  ATI_6 = "data/projects/Around_The_Island_fDOM/around_the_island/run6/",
  ATI_7 = "data/projects/Around_The_Island_fDOM/around_the_island/run7/",
  ATI_8 = "data/projects/Around_The_Island_fDOM/around_the_island/run8/"
)

# Process individual runs -------------

# SGD
SGD_cabral <- process_aqualog(data_directory = run_dirs[["SGD_cabral"]],
                        run_name = "cabral",
                        sample_key_file = "sgd_run_fdom_cabral.tsv",
                        org_file = "cabral_sample_sheet.csv")

SGD_varari <- process_aqualog(data_directory = run_dirs[["SGD_varari"]],
                        run_name = "varari",
                        sample_key_file = "sgd_run_fdom_varari.tsv")


# RRCoral
RRC_run1 <- process_aqualog(data_directory = run_dirs[["RRC_1"]],
                        run_name = "RRCoral_1",
                        sample_key_file = "RRCoral_run_fdom_Sheets - Run1.tsv",
                        org_file = "RRCoral_1_sample_sheet.csv")

RRC_run2 <- process_aqualog(data_directory = run_dirs[["RRC_2"]],
                        run_name = "RRCoral_2",
                        sample_key_file = "RRCoral_run_fdom_Sheets - Run2.tsv",
                        org_file = "RRCoral_2_sample_sheet.csv")

# Spiffy
spiffy <- process_aqualog(data_directory = run_dirs[["Spiffy"]],
                          run_name = "spiffy",
                          sample_key_file = "Spiffy_run_fdom_sheets - Spiffy.tsv")

spiffyLang <- process_aqualog(data_directory = run_dirs[["SL"]],
                              run_name = "spiffy_lang",
                              sample_key_file = "Spiffy_Lang_run_fdom_sheets - Spiffy_Lang.tsv",
                              org_file = "spiffy_lang_sample_sheet.csv")

# Diel
diel_run1 <- process_aqualog(data_directory = run_dirs[["Diel_1"]],
                             run_name = "diel_run1" ,
                             sample_key_file = "diel_run_fdom_sheets - Run1.tsv",
                             org_file= "diel_run1_sample_sheet.csv")

diel_run2 <- process_aqualog(data_directory = run_dirs[["Diel_2"]],
                             run_name = "diel_run2" ,
                             sample_key_file = "diel_run_fdom_sheets - Run2.tsv",
                             org_file= "diel_run2_sample_sheet.csv")

# around the island
ati_run1 <- process_aqualog(data_directory = run_dirs[["ATI_1"]],
                            run_name = "aroundtheisland_run1" ,
                            sample_key_file = "ati_run1.tsv")

ati_run2 <- process_aqualog(data_directory = run_dirs[["ATI_2"]],
                            run_name = "aroundtheisland_run2" ,
                            sample_key_file = "ati_run2.tsv")

ati_run3 <- process_aqualog(data_directory = run_dirs[["ATI_3"]],
                            run_name = "aroundtheisland_run3" ,
                            sample_key_file = "ati_run3.tsv")
# *bad blanks removed
ati_run4 <- process_aqualog(data_directory = run_dirs[["ATI_4"]],
                            run_name = "aroundtheisland_run4" ,
                            sample_key_file = "ati_run4.tsv",
                            org_file = "aroundtheisland_run4_sample_sheet.csv")

ati_run5 <- process_aqualog(data_directory = run_dirs[["ATI_5"]],
                            run_name = "aroundtheisland_run5" ,
                            sample_key_file = "ati_run5.tsv")
# *bad blanks removed
ati_run6 <- process_aqualog(data_directory = run_dirs[["ATI_6"]],
                            run_name = "aroundtheisland_run6" ,
                            sample_key_file = "ati_run6.tsv",
                            org_file = "aroundtheisland_run6_sample_sheet.csv")

# *bad blanks removed
ati_run7 <- process_aqualog(data_directory = run_dirs[["ATI_7"]],
                            run_name = "aroundtheisland_run7" ,
                            sample_key_file = "ati_run7.tsv",
                            org_file = "aroundtheisland_run7_sample_sheet.csv")

ati_run8 <- process_aqualog(data_directory = run_dirs[["ATI_8"]],
                            run_name = "aroundtheisland_run8" ,
                            sample_key_file = "ati_run8.tsv")


all_runs = list(SGD_cabral = SGD_cabral,
                SGD_varari = SGD_varari,
                RRC_run1 = RRC_run1,
                RRC_run2 = RRC_run2,
                spiffy = spiffy,
                spiffyLang = spiffyLang,
                diel_run1 = diel_run1,
                diel_run2 = diel_run2,
                ati_run1 = ati_run1,
                ati_run1 = ati_run2,
                ati_run3 = ati_run3,
                ati_run4 = ati_run4,
                ati_run5 = ati_run5,
                ati_run6 = ati_run6,
                ati_run7 = ati_run7,
                ati_run8 = ati_run8)

saveRDS(all_runs, "data/compiled/all_runs.rds") 




# compile runs -----------------------


# compile run sheets for matlab script
sample_sheets = lapply(all_runs, function(x) x$sample_sheet)


rundt = rbindlist(sample_sheets)


rundt[ , sample := paste(run_name, sample, sep = "/")]
rundt[ , Samplepath := paste(run_name, SEM_file, sep = "/")]
rundt[ , Blank1 := paste(run_name, blank1, sep = "/")]
rundt[ , Blank2 := paste(run_name, blank2, sep = "/")]
rundt[ , abs := paste(run_name, ABS_file, sep = "/")]

rundt[ ,c("blank1","blank2","ABS_file","SEM_file") := NULL]

# check blank count
length(unique(paste(rundt$run_name, rundt$Blank1)))

# check for duplicate sample names
rundt[duplicated(UniqueID), .(UniqueID, run_name, sample)]

# add blanks
for(run in unique(rundt$run_name)){
 
  blanks = unique(rundt[run_name == run, Blank1])
  n_blanks = length(blanks)
  
   dt = data.table(
           UniqueID = paste0("BLANK", 1:n_blanks),
           run_name = run,
           run_order = NA,
           sample = "Blank",
           group = NA,
           Samplepath = blanks,
           Blank1 = c(tail(blanks,-1), head(blanks, 1)),
           Blank2 = c(tail(blanks,-2), head(blanks, 2)),
           abs = sub("BEM", "ABS", blanks)
           )
   rundt = rbindlist(list(rundt, dt))
 }


fwrite(rundt, "data/compiled/matlab_all_sample_sheets.csv")

# run compilation function
compiled_runs = compile_runs(run_dirs = run_dirs,
                             out_dir = "data/compiled")


#  Plot eems -------------------
all_runs = readRDS("data/compiled/all_runs.rds")

# plot samples

for(a_run in all_runs){
  
  a_run_ps = lapply(1:length(a_run$EEMs), function(i){
    plot_eem(a_run$EEMs[[i]], sample_name = names(a_run$EEMs)[i], rows_as_names = T)
    
    })
  
  pdf(paste0("plots/eems/",a_run$indices$run_name[[1]],"_eem_plots.pdf"))
  print(a_run_ps)
  dev.off()
}



# Plot  blanks

for(i in 1:length(run_dirs)){
  
blank = read_eems(file.path(run_dirs[i],"processed_data","blanks"))

bps = lapply(1:length(blank), function(i){
  plot_eem(blank[[i]], sample_name = c(1:20)[i], rows_as_names = F)
  
})

pdf(paste0("plots/eems/blanks/",names(run_dirs)[i],"_blank_plots.pdf"))
print(bps)
dev.off()

}



