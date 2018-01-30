#!/usr/bin/env Rscript
options( warn = -1 )
args = commandArgs(trailingOnly=TRUE)
suppressMessages(source("SPURF.R"))
suppressMessages(require(ggplot2))
suppressMessages(require(ggseqlogo))
AA_LIST = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-')

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "SPURF_output"
}

input_seq <- args[1]
outbase <- args[2]
mode <- if (is.na(args[3])) "l2" else "jaccard"

pred.prof <- predict.prof(input_seq, mode)

mat <- matrix(nrow=21, ncol=149)
for (i in 1:149) {for (j in 1:21) { mat[(i-1)*21+j] <- pred.prof[(i-1)*21+j] }}
rownames(mat) <- AA_LIST
colnames(mat) <- 1:149


offset <- 0
for (i in 1:149) {
  if (mat[21,(i-offset)] == 1) {
    mat <- mat[,-(i-offset)]
    offset <- offset + 1
  }
}

b <- 27
for (i in 1:5) {
  b_str <- as.character(b+i)
  if (sum(b_str == colnames(mat)) == 1) {
    CDR1_min <- (1:ncol(mat))[b_str == colnames(mat)]
    break
   }
}

b <- 41
for (i in 1:5) {
  b_str <- as.character(b-i)
  if (sum(b_str == colnames(mat)) == 1) {
    CDR1_max <- (1:ncol(mat))[b_str == colnames(mat)]
    break
   }
}

b <- 57
for (i in 1:5) {
  b_str <- as.character(b+i)
  if (sum(b_str == colnames(mat)) == 1) {
    CDR2_min <- (1:ncol(mat))[b_str == colnames(mat)]
    break
   }
}

b <- 78
for (i in 1:5) {
  b_str <- as.character(b-i)
  if (sum(b_str == colnames(mat)) == 1) {
    CDR2_max <- (1:ncol(mat))[b_str == colnames(mat)]
    break
   }
}

b <- 108
for (i in 1:5) {
  b_str <- as.character(b+i)
  if (sum(b_str == colnames(mat)) == 1) {
    CDR3_min <- (1:ncol(mat))[b_str == colnames(mat)]
    break
   }
}

b <- 138
for (i in 1:5) {
  b_str <- as.character(b-i)
  if (sum(b_str == colnames(mat)) == 1) {
    CDR3_max <- (1:ncol(mat))[b_str == colnames(mat)]
    break
   }
}

colnames(mat) <- 1:ncol(mat)

suppressMessages(
logo <- ggplot() +
  geom_logo(mat, method = 'prob', seq_type='aa') +
  annotate('segment', x = (CDR1_min-0.5), xend = (CDR1_max+0.5), y=1.1, yend=1.1, size=2) +
  annotate('text', x=(CDR1_min+CDR1_max)/2, y=1.2, label='CDR1', size = 5) +
  annotate('segment', x = (CDR2_min-0.5), xend = (CDR2_max+0.5), y=1.1, yend=1.1, size=2) +
  annotate('text', x=(CDR2_min+CDR2_max)/2, y=1.2, label='CDR2', size = 5) +
  annotate('segment', x = (CDR3_min-0.5), xend = (CDR3_max+0.5), y=1.1, yend=1.1, size=2) +
  annotate('text', x=(CDR3_min+CDR3_max)/2, y=1.2, label='CDR3', size = 5) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels=c('', colnames(mat), '')) +
  theme_classic() +
  theme(text = element_text(size=20), legend.position = "bottom", panel.background = element_rect(fill = NA))
)

write.table(mat, paste(outbase, '.tab', sep=''), sep='\t', quote=F)
suppressMessages(pdf(paste(outbase, '.pdf', sep=''),  width=75/119*ncol(mat), height=3))
suppressMessages(print(logo))
invisible(dev.off())

cat('Wrote results to:\n')
cat(paste(outbase, '.tab\n', sep=''))
cat(paste(outbase, '.pdf\n', sep=''))
