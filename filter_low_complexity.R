# This script is useful for filtering out low complexity sequences.
# It's intended for Illumina amplicon data that contains reads with repeats.
# These reads may have high quality, but the actual base pairs are repeats (e.g. all Gs)
# Low complexity seqs can disrupt dada2 at the learnErrors step.
library(ShortRead)
library(dada2)

# Script parameters ----------------
# S
setwd("/home/sean/Documents/trouble/LearnErrors_troubleshoot/test/")

tm <- proc.time()

# Directory where seq files are located. Should contain all reads (R1 and R2 if paired). 
# Default = current directory.
seq_files_directory = "./"

# Directory where filtered sequences should be written.
output_directory = "cmp_filtered_fastqs"

# Complexity cutoff score (see seqComplexity in dada2 for more information on scoring.
cutoff = 4

# Logical setting for paired reads. If true, removes corresponding reads from paired file during filtering.
# If set to FALSE, treats all reads as R1 reads.
paired_reads = T

# Read 1 pattern. For paired reads, what string pattern in the filename identifies read 1?
R1_pattern = "_R1_"

# Read 2 pattern. For paired reads, what string pattern in the filename identifies read 2?
# Assumes the rest of the filename matches read 1.
R2_pattern = "_R2_"

# Define Functions ----------------------------

# detect_low_comp() determines whether a sequence file has low complexity reads and prints percentage
# specify the complexity score cutoff (5 should be good) and percentage allowable (zero by defautl)
# identify your problem samples!

detect_low_comp <-function(seq_file){
  # get complexity scores
  cmp_scores <-seqComplexity(getSequences(seq_file))
  # determine percentage of reads that are low complexity
  perc_low_comp <- length(cmp_scores[cmp_scores < cutoff]) /
                    length(cmp_scores) 
                   
  
  # if too many low complexity seqs, return file name
  # otherwise, return nothing
  if( perc_low_comp > 0 ){
    print(paste0(seq_file,
                 " Contains low complexity seqs: ",
                 round(perc_low_comp * 100, 2), "%"))
    return(cmp_scores)
    
  }else{
    print(paste0(seq_file,
                 " OKAY"))
    return(NULL)
  }
}

# filter_low_comp() will remove low complexity reads from a fastq file.
# Specify the cutoff for low complexity (11 by default).
# If you are working with paired reads, make sure to specify the other read.
# It's important that the R1 and R2 reads contain the same pairs of sequences.
# R1 can be forward or reverse (e.g. if only filtering reverse reads, set them as R1)
filter_low_comp <- function(R1,
                            paired = paired_reads,
                            cmp_scores = checked_files,
                            outdir = output_directory){
  
  ## Read 1
  # read in the read 1 sequencing file
  raw_R1 <- readFastq(R1)
  # get read 1 complexity scores
  cmps <- cmp_scores[[R1]]
  # generate summary report
  df_out <-  data.frame(Read1 = R1,
                        Total_seqs = length(cmps),
                        Read1_low_comp_seqs = length(cmps[cmps < cutoff]),
                        R1_perc_low_comp = length(cmps[cmps < cutoff]) / length(cmps))
  ## Read 2
  R2 = NULL
  # get complexity scores from paired read, if applicable
  if(isTRUE(paired)){
    # generate read 2 filename based on R1 filename
    R2 <- sub(R1_pattern, R2_pattern, R1)
    
    # check that read 2 file exists
    if(!file.exists(R2)) stop(paste0("ERROR: Read 2 file not found - ", R2))
    
    raw_R2 <- readFastq(R2)
    
    # check read 2 is same length as read 1
    if(length(raw_R2) != length(raw_R1)) stop("ERROR: Read1 and Read2 have different lengths")
    
    # keep whichever complexity score is lower from read 1 and read 2
    R2_cmps <- cmp_scores[[R2]]
    cmps <- mapply(function(x,y) min(x,y), x= cmps, y = R2_cmps)
    
    # add read 2 information to summary report
    df_out$Read2 = R2
    df_out$Read2_low_comp_seqs = length(R2_cmps[R2_cmps < cutoff])
    df_out$Read2_perc_low_comp = length(R2_cmps[R2_cmps < cutoff]) / length(R2_cmps)

  }
  ## Filter Read 1
  filtered_R1 <- raw_R1[cmps > cutoff]
  
  # write out the filtered sequence file
  writeFastq(filtered_R1, file = file.path(outdir, R1))
  
  ## Filter Read 2
  if(isTRUE(paired)){
    ## Filter Read 1
    filtered_R2 <- raw_R2[cmps > cutoff]
    
    # write out the filtered sequence file
    writeFastq(filtered_R2, file = file.path(outdir, R2))
    
    # return summary data.frame
  }
  
  # Print number of seqs filtered
 print(paste0(
       "Filtered low comp from ",
       R1," ", R2
       ))
 print(paste0(
       "Removed: ", length(raw_R1) - length(filtered_R1),
       "  Kept: ", length(filtered_R1) 
 ))
 # return summary report
 return(df_out)
}

# testing
# read1 <- "MV02_S08_S497_L002_R1_trimmed.fastq.gz"
# read2 <- "MV02_S08_S497_L002_R2_trimmed.fastq.gz"
# r1_score <- detect_low_comp(read1)
# r2_score <- detect_low_comp(read2)
# filter_low_comp(read2, read1, cmp_scores = r2_score)


# Run functions on all files in the directory ------------------------

# get a list of files in directory
seq_files <- list.files(path = seq_files_directory)

# make sure the list only includes .fastq or .fq files 
seq_files <- seq_files[grepl("\\.fastq|\\.fq", seq_files)]

# identify read 1 files
R1_files <- seq_files[grepl(R1_pattern, seq_files)]

if(length(R1_files) != length(seq_files)/ 2) warning("Unequal number of R1 and R2 seqs!")
if(length(R1_files) == 0) stop("ERROR: No R1 Seqs! Check pattern and working directory!")

# Get complexity scores for all files -----------
checked_files <- lapply(seq_files,
                         FUN = detect_low_comp)
names(checked_files) <- seq_files

# write out in case next step fails
saveRDS(checked_files, "checked_files.RDS")
# Create output directory for filtered files
if(!exists(output_directory)){
  dir.create(output_directory)
}else{
  warning("Output directory already exists! Remove it before trying again.")
  stop()
}


# Remove low complexity reads from read 1 and, optionally, read 2 -------------
# Call function filter_low_comp() which performs filtering and writes out file to output directory
# The function runs on each read or read pair in the directory.

summary_list <- lapply(R1_files,
                         FUN = filter_low_comp,
                         paired = paired_reads)
# Write out filtering summary
summary_report <- do.call("rbind", summary_list)

write.csv(summary_report, "filter_low_comp_report.csv", row.names = F)

proc.time() - tm
