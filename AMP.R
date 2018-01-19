
predict.prof = function(input_seq) {
  ############ output annotation file path
  output_path = "tmp_annotations.csv"
  
  ## run partis annotation on input sequence
  partis_cmd <- paste('python annotate_sequence.py --sequence', as.character(input_seq),
                      '--outfile', as.character(output_path), sep=" ")
  system(partis_cmd)
  
  df <- read.csv(output_path, header = T)
  input_prof <- df[1, 5:ncol(df)]
  naiveAA <- df[2, 5:ncol(df)]
  neutral_subs <- df[3, 5:ncol(df)]/df$Nseqs[3]
  v_gene = as.character(df$v_gene[1])
  unlink(output_path)
  
  # compile model fitting code
  Rcpp::sourceCpp("stepwise_fitting.cpp")
  
  # load model alphas
  load("cached_data/final_l2_alphas.Rdata")
  alpha = c()
  for (i in 1:4) alpha = c(alpha, alphas[i,])
  
  # load V subgroup/gene profiles
  load("cached_data/public500_vIMGT.Rdata")
  
  # extract the V subgroup/gene identities
  v.subgrp.id = gsub("^IGH(((V[0-9])[/-][NLOR0-9]+([A-Z]|-[0-9]+)?)\\*[0-9]+)(\\+([ACGT][0-9]+[ACGT](\\.[ACGT][0-9]+[ACGT])*|[0-9]+))?$", "\\3", v_gene)
  v.gene.id = gsub("^IGH(((V[0-9])[/-][NLOR0-9]+([A-Z]|-[0-9]+)?)\\*[0-9]+)(\\+([ACGT][0-9]+[ACGT](\\.[ACGT][0-9]+[ACGT])*|[0-9]+))?$", "\\2", v_gene)
  
  # extract the model profiles
  naiveAA.prof = matrix(unlist(naiveAA), nrow=1)
  vgene.prof = matrix(unlist(dff.vgene.centers[v.gene.id,]), nrow=1)
  neut.prof = matrix(unlist(neutral_subs), nrow=1)
  vsubgrp.prof = matrix(unlist(dff.vsubgrp.centers[v.subgrp.id,]), nrow=1)
  input.prof = matrix(unlist(input_prof), nrow=1)
  
  # compute the no-gap profile
  nogap.prof = matrix(NA, nrow=1, ncol=ncol(vsubgrp.prof))
  
  for (j in 1:149) {
    if (all(naiveAA[1, (21*(j-1)+1):(21*j)] == c(rep(0,20), 1))) {
      nogap.prof[1, (21*(j-1)+1):(21*j)] = 0
    } else {
      nogap.prof[1, (21*(j-1)+1):(21*j)] = 1
    }
  }
  
  return(predict_nogap_profs(alpha, rep(1:149, each=21), input.prof,
                             list(naiveAA.prof, vgene.prof, neut.prof, vsubgrp.prof),
                             nogap.prof, 149))
}
