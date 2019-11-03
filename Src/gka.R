gka <- function(response = NA, epsilon = NA){
  # Input: response = vector of responses, 
  #        epsilon = tolerance
  # Output: A summary of GK Algorithm tuples
  
  # Initial setting
  summary <- data.frame(matrix(NA, ncol=3, nrow=1))
  colnames(summary) <- c('v','g','delta')
  s <- 0
  n <- 0
  
  
  # COMPRESS phase
  for (n in 0:(length(response)-1)){
    
    v <- response[n+1]
    
    # Fill in the first 2 iterations
    if (n == 0){
      summary[1,] <- c(v,1,0)
      # Update size of summary
      s <- nrow(summary)
    }
    
    else{

      if (n %% (1/(2*epsilon)) == 0){

        i = s - 1
        while (i >= 2){
          
          j = i-1
          delta.i <- summary[i,3]
          g.sum <- sum(summary[j:i,2])
          v <- summary[i,1]
          
          while (j >= 2 && ((g.sum + delta.i) < 2*epsilon*n)){
            j <- j - 1
            g.sum <- g.sum + summary[j,2]
          }
          
          # Tune one index up
          j <- j + 1
          
          # DELETE phase
          if (j < i){
            # Merge tuples from j to i
            summary <- summary[-((j+1):i),]
            summary[j,] <- data.frame('v'=v, 'g'=g.sum-summary[(j-1),2], 'delta'=delta.i)
          }
          
          # Continue from the largest integer smaller than j
          i <- j - 1
          # Update size of the summary
          s <- nrow(summary)
        }
      }
      
      
      # INSERT phase
      v.0 <- summary[1,1]
      v.s_1 <- summary[s,1]
      
      # Extreme cases
      tuple.new <- data.frame('v'=NA, 'g'=NA, 'delta'=NA)
      if ( v < v.0 ){
        delta <- 0
        new.position <- 0
        summary <- rbind(tuple.new, summary)
      }
      else if ( v > v.s_1 ){
        delta <- 0
        new.position <- s
        summary <- rbind(summary, tuple.new)
      }
      else{
        # Find appropriate index i
        new.position <- which( v < summary[,1] )[1] - 1
        delta <- summary[new.position,2] + summary[new.position,3] - 1
        summary <- rbind(summary, tuple.new)
        summary[(new.position+2):(s+1), ] <- summary[(new.position+1):s, ]
      }
      
      # Insert new tuple
      tuple.new <- data.frame('v'=v, 'g'=1, 'delta'=delta)
      summary[(new.position+1),] <- tuple.new
      # Update size of summary
      s <- s + 1
    }
  }
  
  return(summary)
}