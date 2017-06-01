
#==================================
# TASK 1. Optional step -> remove all data frames in RStudio. Uncomment to use
#==================================  


#rm(list=ls())


#==================================#==================================#===================================
#==================================#==================================#===================================
######### # TASK 2.   			USER SETs  Global Variables: logP-VALUE CUTOFF			##################
#==================================#==================================#===================================
#==================================#==================================#===================================

																						
# set working directory 
#setwd ("C:/Users/etc")						#

#setwd ("E:/Dropbox/Trianta/akadimaika/karolinska_work/DAMATO/manhattan_plot_project/Knitr_shiny_script/knitr")	
#setwd ("E:/Dropbox/Trianta/akadimaika/karolinska_work/DAMATO/manhattan_plot_project/Knitr_shiny_script/knitr/NewDataset/Ek_GWAS_data")
setwd ("C:/Users/tripap/Downloads/newdataset/Ek_GWAS_data")
#Graph settings
##Lines																						
suggestive_cutoff <- 5e-5 # cut off for suggestive line
annotation_cutoff <- 2.3 # log scale.  Set logp cutoff to select those SNPs that get annotated 

##Number of reported SNPs with highest -logP-values. !!! -> KEEP THE "L" after the number (means integer) <- !!!
number_of_SNPs_to_show <- 1000L


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! you need to change the column names to match the dataset you use
# !!! 
# !!! Positions in Code
# !!! 	"colsToKeep" , before dataset initial import, (Task 1. step 4)
# !!! 	"setColumnOrder" , (Task 1. step 3)  and   									!!
# !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################################################







#==================================
# TASK 3. INSTALL/LOAD Prerequisite libraries (only first time) or load directly
#==================================  

if (!require("ggplot2")) install.packages("ggplot2", verbose=FALSE)
if (!require("ggiraph")) install.packages("ggiraph", verbose=FALSE)
if (!require("png")) install.packages("png", verbose=FALSE)
if (!require("grid")) install.packages("grid", verbose=FALSE)
if (!require("gridExtra")) install.packages("gridExtra", verbose=FALSE) # makes nice printed tables
if (!require("sqldf")) install.packages("sqldf", verbose=FALSE)
if (!require("printr")) install.packages("printr", verbose=FALSE)
if (!require("data.table")) install.packages("data.table", verbose=FALSE)
if (!require("dplyr")) install.packages("dplyr", verbose=FALSE)





#==================================
# TASK 4. File import and set working directory
#==================================  
#Step.1 PRE-select columns to keep
#Step.2 Read in using the data
#Step.3 Rearrange columns
#Step.4 Filter-out rows where p=1 to avoid "inf" logp
#==================================  




#Step.1 PRE-select columns to keep
#based on the dataset column names you have	

### !!!! CHANGE COLUMN NAMES DEPENDING ON DATASET!!! <--
#colsToKeep <- c("SNP.UKB_halfcontrols", "CHR", "BP", "OR.meta", "P.meta",  "Gene",  "Dist")

colsToKeep <- c("SNP", "CHR", "BP", "OR", "P")




#Step.2 Read in using the data
# Read in using the data.table method

sink("freadoutput.txt") # Sink / hide output of fread, which shows up even if it should not in knitr


#dat_smp <- fread("dat.smp.txt" , header="auto", na.strings="NA", stringsAsFactors=FALSE, select=colsToKeep, verbose=FALSE)

dat_smp <- fread("daner_IBS_all_0812a_2_copy.txt" , header="auto", na.strings="NA", stringsAsFactors=FALSE, select=colsToKeep, verbose=FALSE)

sink () # return output to the terminal. It is useful if file is parsed in knitr
 
 
# str(dat_smp)






#Step.3 Rearrange columns
### !!!! CHANGE COLUMN NAMES DEPENDING ON DATASET!!! <--
setColumnOrder(dat_smp, c('CHR', 'BP', 'P', 'OR', 'SNP'))

 dat_smp[,3] # Check if P values are on the 3rd column



Buffer_DATA_FRAME <- as.data.frame(dat_smp)  # use the RAW data set as the input. process it as a data.frame
# str(Buffer_DATA_FRAME)





#Step.4 Cutoff rows where p=1 to avoid "inf" logp
Buffer_DATA_FRAME <- Buffer_DATA_FRAME[!(Buffer_DATA_FRAME$P== "1"),]






#==================================
# TASK 5. Calculcations
#==================================
#Step.1 Create duplicate column of original p-vals
#Step.2 Calculate =log10P-val
#Step.3 rename columns
#Step.4 rearrange columns
#==================================  


#Step.1 Create duplicate column of original p-vals
Buffer_DATA_FRAME$pval_original <- Buffer_DATA_FRAME[,3] 						# Keep original p values, put them in a duplicate column "pval_original"


#Step.2 Calculate =log10P-val
Buffer_DATA_FRAME[,3] <- -log10(Buffer_DATA_FRAME[,3]) 				# calc of logP and put these in place of the original P values
Buffer_DATA_FRAME[,5] <- as.character(Buffer_DATA_FRAME[,5]) 		# Set SNP column as characters type
#str(Buffer_DATA_FRAME)



#Step.3 rename columns
colnames(Buffer_DATA_FRAME) <- c('chr', 'pos', 'logp', 'or','SNP', 'original_P')		# rename columns


#Step.4 rearrange columns
Buffer_DATA_FRAME <- Buffer_DATA_FRAME[,c('SNP','chr','pos','logp', 'or','original_P')] 		# rearrange columns

#str(Buffer_DATA_FRAME)



#==================================
# TASK 6.Positioning of ticks in middle of a chromosome, as label of x-axis
#==================================  

#length(unique(Buffer_DATA_FRAME$chr))

chrNum=length(unique(Buffer_DATA_FRAME$chr))

for (i in 1:chrNum){ ndx <- which(Buffer_DATA_FRAME[, 2]==i)  				# for each chr (assingn subset <-  select subeset for each chr)
lstMrk <- max(Buffer_DATA_FRAME[ndx, 3])  									# for this chr. subset, find the last (max) bp position
if (i < chrNum) ndx2 <- which(Buffer_DATA_FRAME[, 2]==i+1)					#if i less than max chr, assign subset for next chr <- select next chr subset

if (i < chrNum) Buffer_DATA_FRAME[ndx2, 3] <- as.numeric(Buffer_DATA_FRAME[ndx2, 3] + lstMrk)	#if i less than max Chr, assign  <- subset of next chr and all BP + add the MAX BP of the previous subset.  The "as.numeric()" has been used to avoid exceeding the limitations of integers of R. without that R would return "NAs produced by integer overflow"


}




bpMidVec <- vector(length=chrNum)

for (i in 1:chrNum){ndx <- which(Buffer_DATA_FRAME[, 2]==i)
posSub <- Buffer_DATA_FRAME[ndx, 3]
bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
}





#==================================
# TASK 7.Prepare graph annotation
#==================================


# prepare the info we want to show when mouse over
Buffer_DATA_FRAME$tips <- paste0("SNP: ",Buffer_DATA_FRAME$SNP, "\n"
                         #,"CHR = ",Buffer_DATA_FRAME$chr, "\n"
                         #,"BP = ",Buffer_DATA_FRAME$pos, "\n"
                         ,"P = ",signif(10^-Buffer_DATA_FRAME$logp,2), "\n"
                         ,"OR = ", round(Buffer_DATA_FRAME$or, 2),sep = " ")
						 
						 



						 
						 
						 
#==================================
# TASK 8.Filter step
#==================================
#Rearrange and sort by logp (highest first-UP).
#Make sure that the top rows have a logP cutoff
#Keep only ____ rows with highest log P values. Number of rows has been set in user settings "number_of_SNPs_to_show"
#
#==================================  



#Rearrange and sort by logp . Highest logp = lowest p val. the (-) REVERSES the order, thus shows from highest to lowest
 order_by_p_val <- Buffer_DATA_FRAME[order(-Buffer_DATA_FRAME$logp),] 

#######################################
# Filter 1.
#######################################

# first find all those above cutoff, then keep the first 1000 Rows
Buffer_DATA_FRAME_head <- ((order_by_p_val[ which(order_by_p_val$logp > annotation_cutoff),])[1:number_of_SNPs_to_show,])



#######################################
# Filter 2. Filter out the results of the "filter 1" from the full dataset
#######################################

# Use results of filter 1, to exclude data (!) from the full data frame 
Buffer_DATA_FRAME_rest <- order_by_p_val[!(order_by_p_val$SNP %in% Buffer_DATA_FRAME_head$SNP),]








					 
#==================================
# TASK 9. Prepare Graph
#==================================

						 
# set some parameters
maxy <- ceiling(max(Buffer_DATA_FRAME$logp)) + 1

# Color palette

pals <- rep(c("#7878ba", "#000041","#7878ba", "#000041","#7878ba", "#000041","#7878ba", "#000041"), 3)[1:length(unique(Buffer_DATA_FRAME$chr))]



#To update default ggplot geom parameters:

#update_geom_defaults("point", list(shape = 16))





#==================================
# TASK 9.1 .Prepare Graph - Build low part, no annotations
#==================================


#build low part with no annotations, based on number of top ____ SNPS

# to control the size of the dots, change the part after the "geom_point"
#p <- ggplot(Buffer_DATA_FRAME_rest, aes(x=pos, y=logp, stroke = 0,  colour=as.factor(chr))) + geom_point(shape = 16 ,  size=1)


p <- ggplot(Buffer_DATA_FRAME_rest, aes(x=pos, y=logp,  colour=as.factor(chr), tooltip = tips), alpha=0.8) + geom_point(shape = 19, size=0.5)

p3 <- p + scale_y_continuous(limits = c(0, maxy),breaks = seq(0, maxy, by = 1)) + scale_color_manual(values=pals)
p4 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() )
p6 <- p4 + scale_x_continuous(labels=as.character(1:chrNum), breaks=bpMidVec) + theme(axis.text=element_text(size=18))
p10 <- p6 + geom_point() + theme_void() + theme(legend.position='none')
ggsave(filename = "manhattan.png", plot = p10, width = 16, height = 9, dpi = 300)



#==================================
# TASK 9.2 .Prepare Graph - Build Upper part, with annotations
#==================================

#build high part with the annotations, based on number of top ____ SNPS

# to control the size of the dots, change the part after the "geom_point"
p <- ggplot(Buffer_DATA_FRAME_head, aes(x=pos, y=logp,   colour=as.factor(chr), tooltip = tips), alpha=0.8) + geom_point(shape = 16, size=0.1)

p3 <- p + scale_y_continuous(limits = c(0, maxy),breaks = seq(0, maxy, by = 1)) + scale_color_manual(values=pals)
p4 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() )
p5 <- p4 + theme(legend.position='none')
p6 <- p5 + scale_x_continuous(labels=as.character(1:chrNum), breaks=bpMidVec) + theme(axis.text=element_text(size=18))
p8 <- p6 + labs(y=expression('-log'[10]*'(p)'), x='Chromosome') + theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))


img <- readPNG("manhattan.png")
p20 <-  p8 + annotation_custom(rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc")), -Inf, Inf, -Inf, Inf) + geom_point_interactive()
p21 <- p20 + geom_point_interactive() + geom_hline(yintercept=-log(suggestive_cutoff, base = 10), linetype=1, col='blue', lwd=1, alpha=1/3) + geom_hline(yintercept=-log(5e-8, base = 10), linetype=1, col='red', lwd=1, alpha=1/3)


ggiraph(code = {print(p21)}, width_svg = 16, height_svg = 9, tooltip_opacity = .90)
# to save as html I just export it as web page in Rstudio...


