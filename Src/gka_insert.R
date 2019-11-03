gka.insert <- function(v=NA,gka.summary=NA){
  # Input: v = new observation
  #        gka.summary = a data frame returned by the gka function
  # Output: the new gka.summary dataframe
  
  # INSERT phase
  s <- nrow(gka.summary)
  v.0 <- gka.summary[1,1]
  v.s_1 <- gka.summary[s,1]
  
  # Extreme cases
  tuple.new <- data.frame('v'=NA, 'g'=NA, 'delta'=NA)
  if ( v < v.0 ){
    delta <- 0
    new.position <- 0
    gka.summary <- rbind(tuple.new, gka.summary)
  }
  else if ( v > v.s_1 ){
    delta <- 0
    new.position <- s
    gka.summary <- rbind(gka.summary, tuple.new)
  }
  else{
    # Find appropriate index i
    new.position <- which( v < gka.summary[,1] )[1] - 1
    delta <- gka.summary[new.position,2] + gka.summary[new.position,3] - 1
    gka.summary <- rbind(gka.summary, tuple.new)
    gka.summary[(new.position+2):(s+1), ] <- gka.summary[(new.position+1):s, ]
  }
  
  # Insert new tuple
  tuple.new <- data.frame('v'=v, 'g'=1, 'delta'=delta)
  gka.summary[(new.position+1),] <- tuple.new
  
  return(gka.summary)

}