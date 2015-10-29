#' Corrige fr. esperadas muy pequeñas
#' @param E, es una lista (de longitud nitems) con fr. esperadas de respuesta incorrectas y correctas para cada item
#' @param O, es una lista (de longitud nitems) con fr. observadas de respuesta incorrectas y correctas para cada item
#' @param mincell, es el minimo numero esperado de personas con score k que contestan o no al item j


collapseCells <- function(O, E, mincell = 1){
  for(i in 1L:length(O)){ #para i entre 1 y nitems
    On <- O[[i]]
    En <- E[[i]]
    
    L <- En < mincell & En != 0
    while(any(L, na.rm = TRUE)){
      if(!is.matrix(L)) break
      whc <- min(which(rowSums(L) > 0L)) #selecciona los scores con problemas
      if(whc == 1L){ #soluciona problemas para el primer score
        En[2L,] <- En[2L, ] + En[1L,]
        On[2L,] <- On[2L, ] + On[1L,]
        En <- En[-1L,]; On <- On[-1L,]
      } else if(whc == nrow(En)){#soluciona problemas para el ultimo score
        En[nrow(En)-1L,] <- En[nrow(En)-1L, ] + En[nrow(En),]
        On[nrow(On)-1L,] <- On[nrow(On)-1L, ] + On[nrow(On),]
        En <- En[-nrow(En),]; On <- On[-nrow(On),]
      } else {#soluciona problemas para scores intermedios
        ss <- c(sum(On[whc-1L,]), sum(On[whc+1L,]))
        up <- (min(ss) == ss)[1L]
        pick <- if(up) whc-1L else whc+1L
        En[pick,] <- En[pick, ] + En[whc,]
        On[pick,] <- On[pick, ] + On[whc,]
        En <- En[-whc,]; On <- On[-whc,]
      }
      L <- En < mincell & En != 0
    }
    En[En == 0] <- NA
    E[[i]] <- En
    O[[i]] <- On
  }###aca termina el for
  return(list("O"=O, "E"=E))
}
