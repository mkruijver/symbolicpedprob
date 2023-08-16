pr_obs_given_ibd_state <- function(ibd_state, observations){

  if (!is.integer(ibd_state)){
    stop("ibd_state needs to be integer")
  }

  number_of_placeholders <- max(ibd_state)

  placeholders <- rep("-1", number_of_placeholders)

  observations
  pr_obs_given_ibd_state_recursive(i_observation = 1L, placeholders = placeholders,
                                   ibd_state, observations)
}

pr_obs_given_ibd_state_recursive <- function(i_observation, placeholders,
                                             ibd_state, observations){

  if (i_observation == (1 + nrow(observations))){
    return("1")
  }
  if (i_observation > nrow(observations)) stop("recursion out of bounds")

  obs_a <- observations[i_observation, 1]
  obs_b <- observations[i_observation, 2]

  i <- as.integer(ibd_state[2*i_observation - 1])
  j <- as.integer(ibd_state[2*i_observation])

  i_free <- placeholders[i] == "-1"
  j_free <- placeholders[j] == "-1"

  # homozygous
  if (obs_a==obs_b){

    i_valid <- i_free || (placeholders[i]==obs_a)
    j_valid <- j_free || (placeholders[j]==obs_a)

    valid <- i_valid && j_valid

    if (valid){

      pr <- "1"

      if (i_free){
        placeholders[i] <- obs_a
        pr <- obs_a
      }
      if (j_free){
        placeholders[j] <- obs_a
        pr <- paste0(pr, "*", obs_a)
      }

      return(paste0(pr, "*",
                    pr_obs_given_ibd_state_recursive(i_observation + 1,
                                                     placeholders,
                                                     ibd_state,
                                                     observations)
      ))
    }
    else{
      return("0")
    }
  }else{

    #heterozygous
    pr_het <- "("

    ## order 1
    i_valid <- i_free || (placeholders[i]==obs_a)
    j_valid <- j_free || (placeholders[j]==obs_b)

    pr_order1 <- "1"

    if (i_valid && j_valid){
      if (i_free){
        # assign i
        placeholders[i] <- obs_a
        pr_order1 <- obs_a
      }

      if (j_free){
        # assign j
        placeholders[j] <- obs_b
        pr_order1 <- paste0(pr_order1, "*", obs_b)
      }

      # recurs
      pr_order1_recurs <- pr_obs_given_ibd_state_recursive(i_observation + 1,
                                                           placeholders,
                                                           ibd_state,
                                                           observations)

      if (pr_order1_recurs == "0"){
        pr_order1 <- "0"
      }
      else if (pr_order1_recurs!="1"){
        pr_order1 <- paste0(pr_order1, "*",
                            pr_order1_recurs)
      }

      # restore state
      if (i_free){
        placeholders[i] <- "-1"
      }
      if (j_free){
        placeholders[j] <- "-1"
      }

    }
    else{
      pr_order1 <- "0"
    }

    ## order 2
    i_valid <- i_free || (placeholders[i]==obs_b)
    j_valid <- j_free || (placeholders[j]==obs_a)

    pr_order2 <- "1"

    if (i_valid && j_valid){
      if (i_free){
        # assign i
        placeholders[i] <- obs_b
        pr_order2 <- obs_b
      }

      if (j_free){
        # assign j
        placeholders[j] <- obs_a
        pr_order2 <- paste0(pr_order2, "*", obs_a)
      }

      # recurs
      pr_order2_recurs <- pr_obs_given_ibd_state_recursive(i_observation + 1,
                                                           placeholders,
                                                           ibd_state,
                                                           observations)

      if (pr_order2_recurs=="0"){
        pr_order2 <- "0"
      }
      else if (pr_order2_recurs!="1"){
        pr_order2 <- paste0(pr_order2, "*",
                            pr_order2_recurs)
      }

      # restore state
      if (i_free){
        placeholders[i] <- "-1"
      }
      if (j_free){
        placeholders[j] <- "-1"
      }

    }
    else{
      pr_order2 <- "0"
    }

    if (pr_order1 == "0" && pr_order2 == "0"){
      pr_het <- "0"
    }
    else if (pr_order1 == "0"){
      pr_het <- pr_order2
    }
    else if (pr_order2 == "0"){
      pr_het <- pr_order1
    }
    else{
      pr_het <- paste0("(", pr_order1, " + ", pr_order2, ")")
    }

    return(pr_het)
  }

}
