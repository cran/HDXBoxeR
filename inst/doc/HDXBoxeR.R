## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center"
)

## ----eval=FALSE---------------------------------------------------------------
#  ####installation using CRAN
#  install.packages("HDXBoxeR") #execute only once
#  

## ----eval=FALSE---------------------------------------------------------------
#  ####installation using Github.
#  #run only if devtools package is not installed on your machine
#  install.packages("devtools")
#  
#  library(devtools) #run next two commends only once
#  devtools::install_github("mkajano/HDXBoxeR")
#  
#  #Once installed, you can load the HDXBoxeR package using the following command:
#  library(HDXBoxeR)
#  

## -----------------------------------------------------------------------------
# Path to example input
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
### Later change the path to the path you your data!!!

##########Deuteration uptake vs percent deuteration file preparation
#### Input for uptake deuteration. Default will use all common states, all common Deuteration Times, max. common number of replicates, no csv will be written, and peptide sequences will be matched. 
input_uptake_timepoint<-output_tp(filepath = file_nm) 
# input for percent deuteration
input_deut_timepoint<-output_tp(filepath = file_nm, percent=T) 

#################
####Time courses 
#Deuteration uptake
input_uptake_timecourse<-output_tc(file_nm) 
# Time courses comparisons
input_deut_timecourse<-output_tc(file_nm, percent=T)

#################
####Analyze selected protein.states, deuteration times, replicates. 
names_states<- nm_states(file_nm) ## returns names of the states in file.
input_two_states_two_timepoints<-output_tp(filepath=file_nm, replicates=3, states=names_states[c(1,2)], times=c("3.00s", "72000.00s"),  percent=FALSE)
###add option to match sequence and to save as csv (here not used)
input_subset_states<-output_tp(filepath=file_nm, replicates=3, states=names_states[2], times=c("3.00s", "72000.00s"), percent=FALSE, seq_match=FALSE, csv="NA")


## -----------------------------------------------------------------------------
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
input_uptake_timepoint<-output_tp(filepath = file_nm) 



###average calculation for all peptides
av1<-ave_timepoint(input_uptake_timepoint)
head(av1)
###average differences (against first Protein State in the file)
da1<-dif_ave(av1)

###standard deviation calculation for all peptides
sd1<-sd_timepoint(input_uptake_timepoint)
head(sd1)

### Global critical interval for 1 protein sets
CI_single(s1 =sd1[,7], replicates = 3 )
##Global critical interval for 2 protein sets
CI_2pts(s1 = sd1[,7], s2=sd1[,8], replicates = 3)

##individual peptide p-value calculation against first set in the input file
pv1<-pv_timepoint(input_uptake_timepoint)
head(pv1)


####Analyze selected protein.states, deuteration times, replicates. 
names_states<- nm_states(file_nm) ## returns names of the states in file.
### choosing different state as a control for analysis
input_states_reversed<-output_tp(file_nm, states=names_states[c(2,1)] )
pv2<-pv_timepoint(input_states_reversed)
head(pv2)

## -----------------------------------------------------------------------------
###load data

states<-arguments_call1(filepath=file_nm)
times<-arguments_call2(filepath=file_nm, states=states)
replicates<-arguments_call3(filepath=file_nm, states=states, times=times)


## -----------------------------------------------------------------------------


#################
####Time courses 
#Deuteration uptake
input_uptake_timecourse<-output_tc(file_nm) 
# Time courses comparisons

small_df1<-input_uptake_timepoint[,c(1:6, 10:12)]

#if one wants to focus just on the few first peptides in the data it is possible to select a few rows from the input data. 
small_df2<-input_uptake_timepoint[1:3,]

#To have more control of the selecting one can also use the select_indices functions as follows:
##for output_tp function variables that can be used are:
#below one time of 60s was used for peptides that started at residue 50, ended at residue 100 and  at the max 12 residues long. 
#not all parameters need to be used.
inda<-select_indices(input_uptake_timepoint,  times = c("60.00s"),start = 50, end=100, length=12)

#after the indices are selected the input can be subsetted as normal

small_df3<-input_uptake_timepoint[inda,]

##for the output_tc function allowed are the following parameters
timecourse_output<-output_tc(file_nm) 
indb<-select_indices(timecourse_output,  states = "bound",start = 50, end=100, length=12)

head(timecourse_output[indb,1:6])


## ----figuptake, fig.height = 3.5, fig.width = 5-------------------------------
###load data
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
input_deut_timepoint<-output_tc(file_nm, percent=TRUE) 

# define the timepoint using in the experiment. Here the times in seconds
x<-c(3, 60,1800, 72000)
uptake_plots(input_deut_timepoint[1:2,],x)


## ----fig1, fig.height = 5, fig.width = 3.5------------------------------------
###load data
# library(HDXBoxeR)
# file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
# input_uptake_timepoint<-output_tp(file_nm) 


### Returns boxplots for all the time points and all protein states. 
boxplot_tp(input_uptake_timepoint, col= c("gold2", "dodgerblue"))

## ----fig2, fig.height = 5, fig.width = 4--------------------------------------
###load data
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
# ## input for deuteration uptakeinput_uptake_timepoint<-output_tp(file_nm) 
# # input for procent deuteration
input_deut_timepoint<-output_tp(file_nm, percent=T) 

## deuteration uptake for time points
plots_av_tp(input_uptake_timepoint)
### average plots with colors chosen by user. 
plots_av_tp(input_uptake_timepoint,replicates=3, cola=c(1:10))
### average percent deuteration
plots_av_tp_proc(input_deut_timepoint)


## ----fig4, fig.height = 5, fig.width = 4--------------------------------------
# ###load data
# file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
# ## input for deuteration uptake
# input_uptake_timepoint<-output_tp(file_nm) 


## difference in uptake deuterations between two sets
plots_diff_tp(input_uptake_timepoint)

#Load deut data
# input for procent deuteration
# input_deut_timepoint<-output_tp(file_nm, percent=T) 

### difference in procent deuteration
plots_diff_tp_proc(input_deut_timepoint,replicates=3, cola=4)

#input with different order states
names_states<- nm_states(file_nm)
input_reversed_states<-output_tp(file_nm, states = rev(names_states))
### average plots for deuteration uptake where control state was chosen differently
plots_diff_tp(input_reversed_states, col="darkgreen")


## ----fig5, fig.height = 5, fig.width = 4--------------------------------------
###load data
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for deuteration uptake
a<-output_tp(file_nm) 

### basic volcano plot
plots_vol_tp(a)
# change colors for significant peptides in volcano plots
plots_vol_tp(a, cola=c(2,3), replicates=3)
#change pv_cutoff from 0.01 to 0.1
plots_vol_tp(a, pv_cutoff = 0.1)


## ----fig6, fig.height = 5, fig.width = 4--------------------------------------
###load data
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for deuteration uptake
a<-output_tp(file_nm) 

### significantly different peptides are colored in red-blue scheme. 

#par(mfrow=c(4,1))  # uncomment if desired
plot_peptide_sig_tp(a,replicates = 3) 

### Plot where 18 peptides per row are drawn (nb_pep_row=18, default=50), 
#p-value&critial interval was made more stingent (0.001)
# % ranges are colored were changed. 
#par(mfrow=c(4,1)) #uncomment if desired
plot_peptide_sig_tp(a,nb_pep_row = 18, ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf), pv_cutoff = 0.001) 


## ----leg1, fig.height = 3, fig.width = 2.5------------------------------------
## Legend for significant peptides plot
#default ranges for figures does not require argument
legend_sig_peptides()
### Using different range scheme
legend_sig_peptides(ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf))


## ----fig7, fig.height = 3, fig.width = 5--------------------------------------
###load data
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for deuteration uptake
a<-output_tp(file_nm) 
# input for procent deuteration
b<-output_tp(file_nm, percent=T) 

##Average uptake heat maps
### heat maps 
plot_heat_map_tp(a, replicates=3, mar_x=3)
##change some parameters
plot_heat_map_tp(a, mar_x=1, 
                 ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf), pv_cutoff=0.01)

### heat map for percent deuteration, require both uptake and percent deuteration data frame as input
plot_heat_map_tp_proc(input_up = a, input_proc = b, replicates=3)

###Maximum uptake or percent deuteration per peptide
plot_heat_map_max_uptake_tp(a, replicates=3)
### #
plot_heat_map_max_uptake_tp_proc(input_up = a, input_proc = b, replicates=3, mar_x=1, ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf), pv_cutoff=0.01)

## ----leg2, fig.height = 3, fig.width = 2.5------------------------------------
## Legend for the heatmaps.
#default ranges for figures does not require argument
legend_heat_map(ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf))


## ----fig_woods, fig.height = 5, fig.width = 5---------------------------------
##
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 


tp_aP<-output_tp(filepath = file_nm, percent=T) #timepoints
indp<-select_indices(tp_aP,  times = c("3.00s", "60.00s"))
deuteration_woods_timepoints(tp_aP[indp,], replicates = 3)

## input for timecourses 
tc_aP<-output_tc(filepath = file_nm, percent=T) #timecourse
deuteration_woods_timecourse(tc_aP, replicates=3)




## ----fig9, fig.height = 4, fig.width = 5--------------------------------------
##
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for timecourses for procent deuteration and uptake
tc_aP<-output_tc(filepath = file_nm, percent=T)

## average procent deuteration in timecourse
plots_av_tcourse(tc_aP, replicates = 3)
###average heat map for timecourses
plot_heat_map_tc(tc_aP, replicates=3)
legend_heat_map_timecourse()



## ----fig8, fig.height = 5, fig.width = 5--------------------------------------
###load data
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for timecourses for procent deuteration and uptake
tc_a<-output_tc(filepath = file_nm)
tc_aP<-output_tc(filepath = file_nm, percent=TRUE)

###All peptides drawn 
robot_plot_All(thP = tc_aP, th=tc_a)
###USe more stingent parameters to have less peptides drawn. CI_factor is a factor that modifies (multiplicate) Critial interval. 
robot_plot_All(thP = tc_aP, th=tc_a, pv_cutoff=0.005, CI_factor = 5)



## ----fig8a, fig.height = 4, fig.width = 5, echo=FALSE-------------------------
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for timecourses for procent deuteration and uptake
tc_a<-output_tc(filepath = file_nm)
tc_aP<-output_tc(filepath = file_nm, percent=T)
names_states<- nm_states(file_nm) ## only two protein states should be used in functions below. 

## returns dataframe with all significant peptides with pvalue=0.01, CI_factor=2
robot_indexes_df(thP = tc_aP, th=tc_a, states=rev(names_states)[1:2], 
                    pvalue = 0.005, CI_factor=1.5)

#list of indexes 
inds<-robot_indexes(thP = tc_aP, th=tc_a, states=rev(names_states)[1:2], 
                    pvalue = 0.005, CI_factor=1.5) ###return indexes of peptide

inds2=inds[c(1,4, 7,8, 9, 11,13,14, 15,16, 20,21, 24)] ### pick which peptides you want to keep 

# Make a final robot plot. Above the plot there is bar that showing lack of coverage in the sets (grey), no coverage on plot (blue), coverage on plot (orange). 
robot_2states_indexes(thP = tc_aP, th=tc_a, states=rev(names_states)[1:2],indexes = inds2, pvalue=0.001, CI_factor = 2, ylim=c(-120, 120), xlim=c(50, 230))


## ----fig_woods_CI, fig.height = 5, fig.width = 5------------------------------
###load data
library(HDXBoxeR)
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 
## input for timecourses for procent deuteration and uptake
tc_a<-output_tc(filepath = file_nm)
tc_aP<-output_tc(filepath = file_nm, percent=TRUE)

###All peptides drawn 
woods_CI_plot(thP = tc_aP, th=tc_a)
###USe more stingent parameters to have less peptides drawn. CI_factor is a factor that modifies (multiplicate) Critial interval. 
woods_CI_plot(thP = tc_aP, th=tc_a, pv_cutoff=0.005, CI_factor = 5)



## ----echo = FALSE, fig10a, fig.height = 4, fig.width = 3, eval=FALSE----------
#  file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#  ## input for deuteration uptake
#  a<-output_tp(file_nm)
#  # input for procent deuteration
#  b<-output_tp(file_nm, percent=T)
#  
#  
#  ### Scripts written for significantly different peptides (uptake data). Color scheme, ranges to be colored and p-value can be changed.
#  #path where the output files will be written needs to be specified!
#  pymol_script_significant_peptide(a, path=tempdir())
#  pymol_script_significant_peptide(a,  ranges=c(-Inf, seq(-120, 120, by=10), Inf), pv_cutoff = 0.005,replicates = 3, order.pep=FALSE, path=tempdir())
#  ###same but for procent deuteration
#  pymol_script_significant_peptide_proc(input_proc = b, input_up = a, path=tempdir(), ranges=c(-Inf, seq(-120, 120, by=10), Inf), pv_cutoff = 0.005,replicates = 3)
#  
#  ###scripts prepared by residue
#  pymol_script_significant_residue(a,  ranges=c(-Inf, seq(-120, 120, by=10), Inf), path=tempdir(), pv_cutoff = 0.005,replicates = 3)
#  

## ----echo = FALSE-------------------------------------------------------------
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR") 

all_summary(file_nm, Dfact=0.85,replicates = 3)


## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#  
#  ###final output for the paper
#  ### USAGE output_prep(pathto_allresults.csv, output_name.csv)
#  output_prep(filepath = file_nm, output_name = tempfile())
#  
#  ###long version of the uptake data, timepoint
#  ### to save file provide a output name to csv flag
#  output_tp(filepath =file_nm , csv = tempfile())
#  
#  
#  ### verbose versions of the output functions create csv files with all important data averages, sd, and pvalues for uptake data
#  verbose_timepoint_output(filepath = file_nm, output_name =tempfile())
#  verbose_timepoint_output(filepath = file_nm, output_name = tempfile(), percent=T)
#  

## ----eval=FALSE---------------------------------------------------------------
#  # returns empty column if data is missing
#  path_to_folders<-system.file("extdata",  package = "HDXBoxeR")
#  
#  extreme_input_gap(hm_dir =path_to_folders, replicates = 3,
#                    timepoints =c(3, 60, 1800, 72000), output_path=tempdir())
#  # if data is missing, it writes there values for undeuterated
#  extreme_input_undeut(hm_dir =path_to_folders, replicates = 2,
#                       timepoints =c(3, 60, 1800, 72000), output_path=tempdir())
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(HDXBoxeR)
#  file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#  
#  ############################
#  ## Load data
#  ###########################
#  
#  ##########Deuteration uptake vs procent deuteration file preparation
#  #### Input for uptake deuteration
#  a<-output_tp(file_nm)
#  # input for percent deuteration
#  b<-output_tp(file_nm, percent=T)
#  
#  # Time courses inputs
#  tc_a<-output_tc(filepath = file_nm)
#  tc_aP<-output_tc(filepath = file_nm, percent=T)
#  
#  
#  ###########################
#  ##definition of some important variables
#  ###################
#  
#  states<-arguments_call1(filepath=file_nm)
#  times<-arguments_call2(filepath=file_nm, states=states)
#  replicates<-arguments_call3(filepath=file_nm, states=states, times=times)
#  
#  ######################################
#  ## Statistical analysis
#  #####################################
#  
#  ###average calculation for all peptides
#  av1<-ave_timepoint(a)
#  ###average differences (against first Protein State in the file)
#  da1<-dif_ave(av1)
#  ###standard deviation calculation for all peptides
#  sd1<-sd_timepoint(a)
#  ### Global critical interval for 1 protein sets
#  CI_single(s1 =sd1[,7], replicates = 3 )
#  ##Global critical interval for 2 protein sets
#  CI_2pts(s1 = sd1[,7], s2=sd1[,8], replicates = 3)
#  ##individual peptide p-value calculation against first set in the input file
#  pv1<-pv_timepoint(a)
#  
#  ## Plots
#  
#  ##uptake plots
#  
#  # define the timepoint using in the experiment. Here the times in seconds
#  x<-c(4, 60,1800, 72000)
#  #par(mfrow=c(1,2))
#  uptake_plots(tc_aP[1:2,],x)
#  
#  
#  # Boxplots
#  boxplot_tp(a, col= c("gold2", "dodgerblue"))
#  
#  # Volcano plots
#  plots_vol_tp(a, cola=c(2,3), replicates=3, pv_cutoff = 0.01)
#  
#  # Significantly different peptides plot across the sequence
#  #par(mfrow=c(4,1)) # uncomment to have 4 panels
#  plot_peptide_sig_tp(a,nb_pep_row = 18, ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf), pv_cutoff = 0.01)
#  
#  ##Average uptake heat maps
#  ### heat maps
#  plot_heat_map_tp(a, mar_x=8, ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf), pv_cutoff=0.01)
#  
#  ### heat map for procent deuteration, require both uptake and procent deuteration data frame as input
#  plot_heat_map_tp_proc(input_up = a, input_proc = b, replicates=3)
#  
#  ###Maximum uptake or procent deuteration per peptide
#  plot_heat_map_max_uptake_tp(a, replicates=3)
#  ###
#  plot_heat_map_max_uptake_tp(input_up = a, input_proc = b, replicates=3, mar_x=8, ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf), pv_cutoff=0.01)
#  
#  ## Legend for the heatmaps.
#  #default ranges for figures does not require argument
#  legend_heat_map(ranges = c(-Inf, -10, -2.5,0,2.5, 10, Inf))
#  
#  ######################
#  ###Timecourses plots
#  ######################
#  
#  ## average procent deuteration in timecourse
#  plots_av_tcourse(tc_aP, replicates = 3)
#  ###average heat map for timecourses
#  plot_heat_map_tc(tc_aP, replicates=3)
#  legend_heat_map_timecourse()
#  
#  ###Robot plots for the timecourses
#  robot_plot_All(thP = tc_aP, th=tc_a, pv_cutoff=0.01, CI_factor = 1)
#  
#  ##############Woods plots
#  ##timecourse and timepoint woods plots
#  ################
#  
#  indp<-select_indices(b,  times = c("3.00s", "60.00s"))
#  deuteration_woods_timepoints(b[indp,], replicates = 3)
#  
#  ## input for timecourses
#  deuteration_woods_timecourse(tc_aP, replicates=3)
#  
#  ###woods plots
#  woods_CI_plot(thP = tc_aP, th=tc_a)
#  ###USe more stingent parameters to have less peptides drawn. CI_factor is a factor that modifies (multiplicate) Critial interval.
#  woods_CI_plot(thP = tc_aP, th=tc_a, pv_cutoff=0.005, CI_factor = 5)
#  
#  
#  #######################
#  ##Pymol scripts
#  #######################
#  
#  pymol_script_significant_peptide(a, path=tempdir(), ranges=c(-Inf, seq(-120, 120, by=10), Inf), pv_cutoff = 0.01,replicates = 3)
#  ###same but for percent deuteration
#  pymol_script_significant_peptide_proc(input_proc = b, input_up = a,path=tempdir(), ranges=c(-Inf, seq(-120, 120, by=10), Inf), pv_cutoff = 0.01,replicates = 3, order.pep=FALSE)
#  
#  ###scripts prepared by residue
#  pymol_script_significant_residue(a, path=tempdir(), ranges=c(-Inf, seq(-120, 120, by=10), Inf), pv_cutoff = 0.01,replicates = 3)
#  
#  ##################
#  ##Summary
#  ##################
#  all_summary(file_nm, Dfact=0.85,replicates = 3)
#  
#  ###############
#  ##Outputs
#  ###############
#  
#  ###final output for the paper
#  ### USAGE output_prep(pathto_allresults.csv, output_name.csv)
#  output_prep(filepath = file_nm, output_name = tempfile())
#  
#  ###long version of the uptake data, timepoint
#  ### USAGE output_tp_csv(pathto_allresults.csv, output_name.csv)
#  output_tp(filepath =file_nm , csv = tempfile())
#  
#  ### verbose versions of the output functions create csv files with all important data averages, sd, and pvalues for uptake data
#  verbose_timepoint_output(filepath = file_nm, output_name = tempfile(), percent=TRUE)
#  
#  ###########################
#  ##Extreme input preparation
#  ###########################
#  
#  # returns empty column if data is missing
#  extreme_input_gap(hm_dir ="filepath", replicates = 3,
#                    timepoints =c(3, 60, 1800, 72000), output_path=tempdir())
#  # if data is missing, it writes there values for undeuterated
#  extreme_input_undeut(hm_dir ="filepath", replicates = 2,
#                       timepoints =c(3, 60, 1800, 72000), output_path=tempdir())
#  
#  

## -----------------------------------------------------------------------------
##All functions with arguments
lsf.str("package:HDXBoxeR")

