#' Symbolic Pedigree Likelihood
#'
#' @param x A `ped` object.
#'
#' @return A `character` with the symbolic likelihood
#' @export
#'
#' @examples
#' # likelihood for a single heterozygous genotype
#' require(pedtools)
#' het <- singleton() |> addMarker(geno="a/b")
#' sLikelihood(het) # 2*a*b
#'
#' # full siblings
#' fs <- nuclearPed(nch = 2) |> addMarker(geno = c(NA, NA, "a/a", "a/a"))
#' sLikelihood(fs) # "0.25*a^4+0.5*a^3+0.25*a^2"
sLikelihood <- function(x){
  if (!isTRUE(length(x$MARKERS) ==1)){
    stop("x needs to have exactly one marker")
  }

  m <- x$MARKERS[[1]]

  allele_labels <- attr(m, "alleles")

  # to compute the likelihood of the observed data,
  # we consider all multi-person IBD states of persons with data
  # and multiply the prior prob with the likelihood given the IBD state

  # determine persons with data
  idx_observed_rows <- which(apply(m,1, function(row) all(row!=0)))
  observed_persons <- x$ID[idx_observed_rows]

  # convert observations from markers to character matrix
  # each row is a genotype for an observed person
  observations <- matrix(allele_labels[m[idx_observed_rows,]], ncol=2,
                         dimnames = list(observed_persons))

  if (length(observed_persons) > 1){
    # obtain multi-person ibd states for observed persons
    multi_person_IBD <- ribd::multiPersonIBD(x, observed_persons)

    # convert multi_person_IBD to integer matrix for easy indexing
    id_columns_to_integer <- function(id) multi_person_IBD[[id]] |>
      strsplit(fixed = TRUE, split = " ") |>
      c(recursive = TRUE) |> as.integer() |> matrix(ncol=2, byrow=TRUE)

    multi_person_IBD_matrix_int <- do.call(cbind,
                                           lapply(observed_persons, id_columns_to_integer))
    colnames(multi_person_IBD_matrix_int) <- paste0(rep(colnames(multi_person_IBD)[-1], each = 2), c("_a","_b"))


    # obtain symbolic likelihood given each IBD state
    multi_person_IBD$ProbObsSymbolic <- apply(multi_person_IBD_matrix_int, 1,
                                              function(ibd_state){
      pr_obs_given_ibd_state(ibd_state, observations)
    })

    # compute likelihood as sum over IBD state of pr(IBD state) * pr(Obs|IBD state)
    multi_person_IBD$LikSummand <- paste0(multi_person_IBD$Prob, "*", multi_person_IBD$ProbObsSymbolic)
    symbolic_lik_unsimplified <- paste0(multi_person_IBD$LikSummand, collapse = " + ")

    symbolic_lik_simplified <- Ryacas::yac_str(paste0("Simplify(", symbolic_lik_unsimplified, ")"))

    return(symbolic_lik_simplified)
  }
  else if (length(observed_persons) == 1){
    a <- observations[1,1]
    b <- observations[1,2]

    inbreeding_coef <- ribd::inbreeding(x, ids = observed_persons)

    if (a==b){
      # homozygous case
      outbred_term <- paste0(a,"^2")

      if (inbreeding_coef == 0){
        return(outbred_term)
      }
      else{
        inbred_term <- a

        return(paste0(inbreeding_coef,"*", inbred_term, " + ", 1.0 - inbreeding_coef, "*", outbred_term))
      }
    }
    else{
      # heterozygous case

      outbred_term <- paste0("2*", a, "*", b)

      if (inbreeding_coef == 0){
        return(outbred_term)
      }
      else{
        return(paste0(1.0 - inbreeding_coef, "*", outbred_term))
      }
    }

  }
  else{
    return("1")
  }
}
