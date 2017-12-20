
### Hardcoded directory where the script lives:
script.dir <- '/fh/fast/matsen_e/kdavidse/AMP'
setwd(script.dir)
############
annotator_path <- paste(script.dir, '/annotate_sequence.py', sep='')
output_path = paste(script.dir, '/tmp_annotations.csv', sep='')



### This should be an input:
input_seq <- 'CGCAGGACTGTTGANGCCTTCGGAGACCCTGTCCCTCACCTGCGTTGTCTCTGGCGGGTCCTTCAGTGATTACTACTGGAGCTGGATCCATCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGGAGCACCAACTACAACCCGTCCCTCGAAAGTCGAGCCACCATATCAGTAGACACGTCCCAGAACAACCTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACTCGGCTGTGTATTACTGTGCGAGAGGCCCGACTACAATGGCTCACGACTTTGACTACTGGGGCCAGGGAACCCTGGTCACC'


callAnnotator <- function(annotator_path, input_seq, output_path) {
    command <- paste('python', annotator_path, '--sequence', as.character(input_seq), '--outfile', as.character(output_path))
    # print(command)
    system(command)
    return(0)
}
return_code <- callAnnotator(annotator_path = annotator_path, input_seq = input_seq, output_path = output_path)

df <- read.csv(output_path, header = T)
input_seq <- df[1, 5:ncol(df)]
naiveAA <- df[2, 5:ncol(df)]
neutral_subs <- df[3, 5:ncol(df)]/df$Nseqs[3]
one_entry <- list('v_gene'=as.character(df$v_gene[1]), 'd_gene'=as.character(df$d_gene[1]), 'j_gene'=as.character(df$j_gene[1]), 'input_seq'=input_seq, 'naiveAA'=naiveAA, 'neutral_subs'=neutral_subs)
print(one_entry)

unlink(output_path)



