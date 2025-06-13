
# PAIR COUNT FOR IRREGULAR LATTICES

pair_count_IR = function(data, adj.mat, missing.cat = NULL)
{
  datavec = c(data)
  split = as.integer(seq(1,length(datavec)+1, length.out = 20))
  
  output.table = data.frame(couple = factor(),
                            abs.frequency = integer(),
                            proportion = numeric())
  
  for(i in 1:(length(split)-1))
  {
    couplevec = NULL
    for(j in split[i]:(split[i+1]-1))
    {
      ind = which(adj.mat[j,]==1)
      
      if (length(ind)>0) 
      {
        assoc = datavec[ind]
        assoc.num = as.numeric(as.factor(datavec))[ind]
        
        couple = ifelse(as.numeric(as.factor(datavec))[j] <= assoc.num,
                        paste(datavec[j], assoc, sep = ""),
                        paste(assoc, datavec[j], sep = ""))} else couple = NULL
                        couplevec = c(couplevec, couple)
    }
    
    if (is.numeric(c(data)))
      cat.names = sort(c(unique(c(data)), missing.cat))
    
    if (is.character(c(data))|is.factor(c(data)))
      cat.names = sort(as.factor(c(unique(c(data)), as.character(missing.cat))))
    couple.names = NULL
    
    for(i in 1:length(cat.names))
      couple.names = c(couple.names, paste(cat.names[i], cat.names[i:length(cat.names)], 
                                           sep = ""))
    couple.n = choose(length(datavec[!is.na(datavec)]),2)
    
    # Build relative frequencies
    couple.list = sort(unique(couplevec))
    abs.freq = as.numeric(table(couplevec))
    abs.freq.complete = numeric(length(couple.names))
    
    for(cc in 1:length(couple.names))
    {
      which.ind = which(couple.list == couple.names[cc])
      if (length(which.ind) > 0) abs.freq.complete[cc] = abs.freq[which.ind]
    }
    rel.freq = abs.freq.complete/couple.n
    freq.table = data.frame(couple.names, abs.freq.complete, rel.freq)
    colnames(freq.table) = c("couple", "abs.frequency", "proportion")
    
    output.table = bind_rows(freq.table, output.table) %>%
      group_by(couple) %>%
      summarise_all(list(~sum(., na.rm = TRUE)))
  }
  output.table$proportion = as.numeric(output.table$abs.frequency/couple.n)
  return(list(probabilities = output.table, Qk = couple.n))
}

################################################################################

# SHANNON Z FOR IRREGULAR LATTICES

shannonZ_IR = function(data, missing.cat = NULL) 
{
  if(!is.matrix(data) & !is.vector(data)) print("Data must be a matrix or a vector")
  
  datavec = c(data)
  
  indx = c()
  for(i in 1:(length(datavec)-1)){
    if(is.na(datavec[i]) & is.na(datavec[i+1]) & is.na(datavec[i+1]==datavec[i]+1)){
      indx = c(indx,i)
    }
  }
  datavec = datavec[-indx]
  
  if(is.na(datavec[1])) datavec = datavec[-1]
  if(is.na(datavec[length(datavec)])) datavec = datavec[-length(datavec)]
  
  print("Building Z variable...")
  
  indx = which(is.na(datavec))
  adj.mat = matrix(0, length(datavec), length(datavec))
  
  for(j in 1:(length(datavec)-1)) adj.mat[j, (j+1):length(datavec)] = 1
  for(j in 1:length(indx)) adj.mat[indx[j], (indx[j]+1):length(datavec)] = 0
  
  output = pair_count_IR(datavec, adj.mat, missing.cat)
  
  print("Computing entropy...")
  prop = output$probabilities$proportion[output$probabilities$proportion > 0]
  localH = sum(prop*log(1/prop))
  
  print("DONE")
  return(list(probabilities = output$probabilities, shannon.Z = localH))
}


##########################################################################

# SPATIAL ENTROPY IRREGULAR LATTICES

spat_entropy_IR = function(data, adj.list, shannZ, missing.cat=NULL)
{
  ##ingredients:
  #1) Z marginal frequencies
  P.zr = shannZ$probabilities$proportion
  names(P.zr) = shannZ$probabilities$couple
  P.zr
  
  #2) W marginal frequencies and Z|wk conditional frequencies
  n.dist = length(adj.list)
  QQ = sum(shannZ$probabilities$abs.frequency)
  P.zr.cond.wk = vector("list", n.dist)
  P.wk = numeric(n.dist)
  datavec = c(data[!is.na(data)])
  for (dd in 1:n.dist)
  {
    output = pair_count(datavec, adj.list[[dd]], missing.cat)
    P.zr.cond.wk[[dd]] = output$probabilities$proportion
    names(P.zr.cond.wk[[dd]]) = output$probabilities$couple
    P.wk[dd] = output$Qk/QQ
  }
  
  ##partial terms
  res.local = mut.local = numeric(n.dist)
  for(dd in 1:n.dist)
  {
    cond.probs = as.numeric(P.zr.cond.wk[[dd]][P.zr.cond.wk[[dd]] > 0])
    marg.probs = as.numeric(P.zr[P.zr.cond.wk[[dd]] > 0])
    res.local[dd] = sum(cond.probs*log(1/cond.probs))
    mut.local[dd] = sum(cond.probs*log(cond.probs/marg.probs))
  }
  res.global = sum(P.wk*res.local)
  mut.global = sum(P.wk*mut.local)
  
  ##output
  return(list(mut.global = mut.global, res.global = res.global, shannZ = shannZ$shannon.Z,
              mut.local = mut.local, res.local = res.local,
              pwk = P.wk, pzr.marg = P.zr, pzr.cond = P.zr.cond.wk,
              Q = QQ, Qk = P.wk*QQ))
}

