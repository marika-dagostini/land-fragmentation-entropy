##################################
#####     FRAGMENTATION     ######
##################################

dataBL = read.csv("Belluno/data_length.csv")

comuni = unique(dataBL$COMUNE[order(dataBL$COMUNE)])[-1]
comuni

dataBL$COMUNE = as.character(dataBL$COMUNE)
dataBL$COMUNE = ifelse(dataBL$COMUNE == "ArsiÃfÂ¨",'Arsiè', dataBL$COMUNE)
dataBL$COMUNE = ifelse(dataBL$COMUNE == "San NicolÃfÂ² di Comelico", 'San Nicolò di Comelico', dataBL$COMUNE)
dataBL$COMUNE = ifelse(dataBL$COMUNE == "ZoppÃfÂ¨ di Cadore" , 'Zoppè di Cadore', dataBL$COMUNE)
dim(dataBL)
head(dataBL)

comuni = unique(dataBL$COMUNE[order(dataBL$COMUNE)])[-1]
print(comuni)

dataBL_c = read.csv("Belluno/dataBL_100C_all.csv")

dataBL_c$COMUNE = as.character(dataBL_c$COMUNE)
dataBL_c$COMUNE = ifelse(dataBL_c$COMUNE == "ArsiÃfÂ¨",'Arsiè', dataBL_c$COMUNE)
dataBL_c$COMUNE = ifelse(dataBL_c$COMUNE == "San NicolÃfÂ² di Comelico", 'San Nicolò di Comelico', dataBL_c$COMUNE)
dataBL_c$COMUNE = ifelse(dataBL_c$COMUNE == "ZoppÃfÂ¨ di Cadore" , 'Zoppè di Cadore', dataBL_c$COMUNE)

dataBL = left_join(dataBL, dataBL_c[,c(1,10,12)], by = c('id', 'COMUNE'))

dataBL$frag = ifelse(dataBL$DNZone <200,1,0)
dataBL$frag = ifelse(dataBL$Ritagliato_build_pc > 1, 1, dataBL$frag)
dataBL$frag = ifelse(dataBL$LENGTH > 1, 1, dataBL$frag)

shp_comuni = readOGR(dsn = 'Belluno', layer = 'data18_BL_100')

shp_comuni@data$COMUNE = ifelse(shp_comuni@data$COMUNE == 'ArsiÃfÂ¨',
                                'Arsiè', shp_comuni@data$COMUNE)
shp_comuni@data$COMUNE = ifelse(shp_comuni@data$COMUNE == 'San NicolÃfÂ² di Comelico', 
                                'San Nicolò di Comelico', shp_comuni@data$COMUNE)
shp_comuni@data$COMUNE = ifelse(shp_comuni@data$COMUNE == 'ZoppÃfÂ¨ di Cadore', 
                                'Zoppè di Cadore', shp_comuni@data$COMUNE)

sort(unique(shp_comuni@data$COMUNE)[1:61])
shp_comuni@data = left_join(shp_comuni@data[,c(1,7)], dataBL, by = c('id'='id','COMUNE'='COMUNE' ))
dim(shp_comuni@data)

rows = 904
cols = 802 # 725008

for(i in 1:length(comuni)){
  
  print(c(i, comuni[i]))
  
  # # Data Matrix
  grid = dataBL
  grid$frag[grid$COMUNE != comuni[i]] = NA
  grid = grid[order(grid[,'id'], grid[,'frag']),]
  grid = grid[!duplicated(grid$id),]
  
  grid = matrix(grid$frag, nrow = rows, ncol = cols, byrow = FALSE)
  grid = remove_empty(grid, which = c("rows", "cols"), quiet = TRUE)
  # 
  # # Shannon Z
  shZ = shannonZ_IR(grid, breaks = 1000)
  # gc()
  write.csv(do.call("rbind", list(unlist(shZ))), 
            file = paste('Belluno/Results/ShannonZ_frag/shZBL_', comuni[i],'.csv', sep = ''))
  
  grid_comune = subset(shp_comuni, COMUNE == comuni[i])
  grid_comune@data = as.data.frame(grid_comune@data[,8])
  names(grid_comune@data) = 'frag'
  
  # jpeg(filename = paste('Belluno/Results/Maps 100m Frag/Belluno.f_', comuni[i],'.100m.jpg', sep = ''),
  #      width = 1000, height = 1000, quality = 100, res = 120)
  # 
  # plot(grid_comune, col = grid_comune@data$frag, main = comuni[i],
  #      sub = 'Pixel lattice with 100x100m resolution of fragmented areas (Black)')
  # 
  # dev.off()
  
  win = as.owin(grid_comune)
  
  # Centroids
  cc = gCentroid(grid_comune, byid = TRUE)
  coord = as.data.frame(cc@coords)
  ccP = ppp(coord$x, coord$y, win)
  
  # Distance Matrix
  dmat = round(pairdist(ccP), digits = 0)
  dmat[lower.tri(dmat, diag = TRUE)] = NA
  
  # Altieri Entropy
  maxdist = sqrt(diff(win$xrange)^2+diff(win$yrange)^2)
  distbreaks = c(0, 100, 201, 401, 601, 1001, 1501, 2001, maxdist)
  
  shZ = read.csv(paste("Belluno/Results/ShannonZ_frag/shZBL_",comuni[i],".csv", sep = ''))
  
  win = as.owin(grid_comune)
  # Centroids
  cc = gCentroid(grid_comune, byid = TRUE)
  coord = as.data.frame(cc@coords)
  
  maxdist = sqrt(diff(win$xrange)^2+diff(win$yrange)^2)
  distbreaks = c(0, 100, 201, 401, 601, 1001, 1501, 2001, maxdist)
  
  datavec = c(grid[!is.na(grid)])
  
  output.table = data.frame(couple.00 = rep(0,8),
                            couple.01 = rep(0,8),
                            couple.11 = rep(0,8))
  
  for (i in 1:nrow(coord)) {
    print(paste(i,'out of',nrow(coord), sep = ' '))
    for (j in 1:nrow(coord)) {
      dist = sqrt((coord$x[i]-coord$x[j])^2+(coord$y[i]-coord$y[j])^2)
      
      for (dd in 1:(length(distbreaks)-1)) {
        if(dist > distbreaks[dd] & dist <= distbreaks[dd+1]){
          if(datavec[i] == 0 & datavec[j] == 0) 
            output.table[dd,1] = output.table[dd,1]+1
          if(datavec[i] == 0 & datavec[j] == 1) 
            output.table[dd,2] = output.table[dd,2]+1
          if(datavec[i] == 1 & datavec[j] == 0) 
            output.table[dd,2] = output.table[dd,2]+1
          if(datavec[i] == 1 & datavec[j] == 1) 
            output.table[dd,3] = output.table[dd,3]+1
        }
      }
      
    }
    
  }
  
  output.table = output.table/2
  
  P.zr = shZ[8:10]
  names(P.zr) = c('00','01','11')
  P.zr
  
  #2) W marginal frequencies and Z|wk conditional frequencies
  n.dist = length(distbreaks)-1
  QQ = sum(shZ[5:7])
  P.zr.cond.wk = vector("list", n.dist)
  P.wk = numeric(n.dist)
  
  for (dd in 1:n.dist) {
    P.zr.cond.wk[[dd]] = output.table[dd,]/sum(output.table[dd,])
    names(P.zr.cond.wk[[dd]]) = c('00','01','11')
    P.wk[dd] = sum(output.table[dd,])/QQ
  }
  
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
  altieri = list(mut.global = mut.global, res.global = res.global, shannZ = shZ$shannon.Z,
                 mut.local = mut.local, res.local = res.local,
                 pwk = P.wk, pzr.marg = P.zr, pzr.cond = P.zr.cond.wk,
                 Q = QQ, Qk = P.wk*QQ)
  
  i=18
  write.csv(do.call("rbind", list(unlist(altieri))),
            file = paste('Belluno/Results/Altieri_sprawl/AltieriBL_', comuni[i],'.csv',
                         sep = ''))
  
  local.sum = altieri$res.local + altieri$mut.local
  
  barplot(height = rbind(altieri$mut.local/local.sum, altieri$res.local/local.sum),
          beside = FALSE,
          col = c('darkgrey', 'white'),
          names.arg = c('w1','w2',' w3', 'w4', 'w5', 'w6', 'w7','w8'),
          main = comuni[i],
          sub = 'Partial Information (Gray) and Partial Residual Entropy (White)
        in proportional terms for each distance range')
  
  gc()
  
  
  jpeg(filename = paste('Belluno/Results/Plot Altieri sprawl/BellunoEntropy_f_',
                        comuni[i],'_100m.jpg', sep = ''),
       width = 1000, height = 700, quality = 100, res = 120)
  
  barplot(height = rbind(altieri$mut.local/local.sum, altieri$res.local/local.sum),
          beside = FALSE,
          col = c('darkgrey', 'white'),
          names.arg = c('w1','w2',' w3', 'w4', 'w5', 'w6', 'w7','w8'),
          main = comuni[i],
          sub = 'Partial Information (Gray) and Partial Residual Entropy (White)
        in proportional terms for each distance range')
  
  dev.off()
  gc()