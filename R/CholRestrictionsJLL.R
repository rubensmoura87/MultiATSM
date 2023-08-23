#' Impose the zero-restrictions on the Cholesky-factorization from JLL-based models.

#'@param SigmaUnres unrestricted variance-covariance matrix (K X K)
#'@param M          number of domestic unspanned factors per country (scalar)
#'@param G          number of global unspanned factors (scalar)
#'@param N         number of country-specific spanned factors (scalar)
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "none"
#'
#'
#'@keywords internal
#'@return   restricted version the Cholesky factorization matrix from JLL-based models (K x K)
#


CholRestrictionsJLL <- function(SigmaUnres, M, G, N, Economies, DomUnit){

  C <- length(Economies)

  if (DomUnit != "None"){
    IdxDomUnit <- which(DomUnit== Economies) # Index of the dominant country
  }

  # Transform the matrix to be the Cholesky form
  SigmaUnres <- t(chol(SigmaUnres))

  # i) Zero restrictions of global variables on spanned factors
  idx0Global <- G + M
  for (h in 1:C){
    idx1Global <- idx0Global+ N
    SigmaUnres[(idx0Global+1):idx1Global, 1:G] <- 0
    idx0Global <- idx1Global + M
  }

  # ii) Zero restrictions of macro domestic variables on spanned factors
  for (i in 1:C){
    idx0RowMacroSpanned <- G + M
    idx0ColMacroSpanned <- G +(i-1)*(M+N)
    idx1ColMacroSpanned <- idx0ColMacroSpanned + M
    for (h in 1:C){ # Fix the columns and loop through the rows
      idx1RowMacroSpanned <- idx0RowMacroSpanned + N
      SigmaUnres[(idx0RowMacroSpanned+1):idx1RowMacroSpanned, (idx0ColMacroSpanned+1):idx1ColMacroSpanned] <- 0
      idx0RowMacroSpanned <- idx1RowMacroSpanned + M
    }
  }


  # iii) Zero restrictions of spanned factors on macro domestic variables
  for (i in 1:C){
    idx0RowSpannedMacro <- G
    idx0ColSpannedMacro <- G + M +(i-1)*(M+N)
    idx1ColSpannedMacro <- idx0ColSpannedMacro + N
    for (h in 1:C){ # Fix the columns and loop through the rows
      idx1RowSpannedMacro <- idx0RowSpannedMacro + M
      SigmaUnres[(idx0RowSpannedMacro+1):idx1RowSpannedMacro, (idx0ColSpannedMacro+1):idx1ColSpannedMacro] <- 0
      idx0RowSpannedMacro <- idx1RowSpannedMacro + N
    }
  }


  # iv) Zero restrictions of Macro country i on Macro country j
  if ( DomUnit != "None"){
    for (i in 1:C){
      if (i!=IdxDomUnit){
        idx0RowMacroMacro <- G
        idx0ColMacroMacro <- G +(i-1)*(M+N)
        idx1ColMacroMacro <- idx0ColMacroMacro + M
        for (h in 1:C){ # Fix the columns and loop through the rows
          idx1RowMacroMacro <- idx0RowMacroMacro + M
          if (i != h){
            SigmaUnres[(idx0RowMacroMacro+1):idx1RowMacroMacro, (idx0ColMacroMacro+1):idx1ColMacroMacro] <- 0
          }
          idx0RowMacroMacro <- idx1RowMacroMacro + N
        }
      }
    }
  }

  # v) Zero restrictions of Spanned factors of country i on Spanned factors country j
  if (DomUnit != "None"){
    for (i in 1:C){
      if ( i!=IdxDomUnit){
        idx0RowSpannedSpanned <- G + M
        idx0ColSpannedSpanned <- G + M +(i-1)*(M+N)
        idx1ColSpannedSpanned <- idx0ColSpannedSpanned + N
        for (h in 1:C){ # Fix the columns and loop through the rows
          idx1RowSpannedSpanned <- idx0RowSpannedSpanned + N
          if (i != h){
            SigmaUnres[(idx0RowSpannedSpanned+1):idx1RowSpannedSpanned, (idx0ColSpannedSpanned+1):idx1ColSpannedSpanned] <- 0
          }
          idx0RowSpannedSpanned <- idx1RowSpannedSpanned + M
        }
      }
    }
  }

  CholJLL<- SigmaUnres

  return(CholJLL)

}
