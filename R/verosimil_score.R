
#' Verosimilitudes de los scores clasicos para calcular fr. esperadas en Orlando
#' Calcula las veorimilitudes
#' @param pr: matriz de probabilidad
#' @param nitems, cantidad de items considerados en la matriz de probabilidad

s_ss=function(pr,nitems){
  sact = matrix(0,ncol = nitems +1,nrow = G)   ###con los cambios en los indices se parece mucho mas a mirt
  for(m in 1:G){                 
    sant = rep(0,(nitems))         
    sant[1] = 1 - pr[m,1]    #(6)
    sant[2] = pr[m,1]       #(7)
    
    for(k in 2:(nitems)){  #item "i" añadido(empieza desde 2 por q ya incluyo el primer item en (6) y (7) )     
      
      sact[m,1] = (1-pr[m,k]) * sant[1]   #(8)
      
      for(kk in 2:k){ #scores 1:(i-1)
        sact[m,kk] = pr[m,k] * sant[kk-1] + (1 - pr[m,k]) * sant[kk]  #(9)
      }
      sact[m,(k+1)] = pr[m,k] * sant[k] #score "i"  #(10)
      sant = sact[m,]
    }
    
    
  }
  return(sact)
}
