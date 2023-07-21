########################################################################
# implementation of ST-DBSCAN for multidimensional data                #
# Changlock Choi                                                       #
# 2021-05-25                                                           #
########################################################################

########################################################################
# INPUTS :                                                             #
# x     = x-axis coordinate of point data                              #                                  
# y     = y-axis coordinate of point data                              #
# time  = time coordinate of point data                                #
# value = additional attributes of point data for clustering           #
# eps = distance maximum for longitude and latitude                    #
# eps2 =  distance maximum for date                                    #
# eps3 = distance maximum for value                                    #
# minpts = number of points to consider a cluster                      #
########################################################################

mdstdbscan <- function (x, 
                        y, 
                        time,
                        value, 
                        eps, 
                        eps2, 
                        eps3, 
                        minpts) { 

  # Construct dataframe from input vectors
  distdata <- cbind.data.frame(x, y)
  time <- time
  value <- value
   # Get number of rows (data points)
  n <- nrow(distdata)
  
  # Initialize vectors for cluster labels, isseed flag, and cluster count
  classn <- cv <- integer(n)
  # is seed flag
  isseed <- logical(n)
  cn <- integer(1)
  
  for (i in 1:n) {
    # Identify unclassified points
    unclass <- (1:n)[cv < 1]
    
    ##making distance
    a <- data.frame(x = distdata[i, 1], y = distdata[i, 2])
    fordist <- cbind.data.frame(a, distdata)
    idist <- abs(sqrt((fordist[,1] - fordist[,3])^2 + (fordist[, 2] - 
                                                         fordist[, 4])^2))
    # calculate value distance
    forvaluedist <- cbind.data.frame(value[i], value)
    ivaluedist <- abs(forvaluedist[, 1] - forvaluedist[, 2])

    # calculate time distance
    fortime <- cbind.data.frame(time[i], time)
    itimedist <- abs(fortime[, 1] - fortime[, 2])

    # if point i is unclassified
    if (cv[i] == 0) {
      # Define neighbors within the distance bounds on each dimension (spatial, temporal, and value)
      reachables <- intersect(unclass[idist[unclass] <= eps],  
                              unclass[itimedist[unclass] <= eps2])
      reachables <- intersect(reachables, unclass[ivaluedist[unclass] <= eps3])

      # mark the point as noise if the number of neighbors is less than the minimum points
      if (length(reachables) + classn[i] < minpts)
        cv[i] <- (-1)                    
      else {
        # else if the point has enough neighbors
        # start a new cluster and mark the point as a seed (core point)
        cn <- cn + 1                   
        cv[i] <- cn
        isseed[i] <- TRUE
        
        # Remove current point from the neighbors
        reachables <- setdiff(reachables, i)
        
        # Remove current point from the list of unclassified points
        unclass <- setdiff(unclass, i) 

        # Increase cluster count for the neighbors
        classn[reachables] <- classn[reachables] + 1
        
        # While there are reachable points (neighbors),
        # the loop keeps expanding the cluster by including reachable points that have enough neighbors
        while (length(reachables)) {

          # Assign current cluster label to reachable points
          cv[reachables] <- cn  
          
          # Points to consider in this iteration
          ap <- reachables  
          
          # Clear reachables for next iteration
          reachables <- integer()

          # for all reachable points j
          for (i2 in seq(along = ap)) {
            j <- ap[i2]
            
            ##make distance again when cluster is expanding
            # Compute distances for point 'j' in the current iteration
            b <- data.frame(x = distdata[j, 1], y = distdata[j, 2])
            jfordist <- cbind.data.frame(b, distdata)
            jdist <- sqrt((jfordist[,1] - jfordist[,3])^2 + 
                            (jfordist[, 2] - jfordist[, 4])^2)
            # value distance
            jforvaluedist <- cbind.data.frame(value[j], value)
            jvaluedist <- abs(jforvaluedist[, 1] - jforvaluedist[, 2])

            # time distance
            jfortime <- cbind.data.frame(time[j], time)
            jtimedist <- abs(jfortime[, 1] - jfortime[, 2])

            # identify neighors for j
            jreachables <- intersect(unclass[jdist[unclass] <= eps],  
                                     unclass[jtimedist[unclass] <= eps2])
            jreachables <- intersect(jreachables, unclass[jvaluedist[unclass]
                                                          <= eps3])
            # If point 'j' has enough neighbors, mark it as a seed and update cluster assignment
            if (length(jreachables) + classn[j] >= minpts) {
              isseed[j] <- TRUE
              cv[jreachables[cv[jreachables] < 0]] <- cn
              reachables <- union(reachables, jreachables[cv[jreachables] == 0])
            }

            # Update cluster size
            classn[jreachables] <- classn[jreachables] + 1

            # Remove the processed point from the list of unclassified points
            unclass <- setdiff(unclass, j)
          }
        }
      }
    }
     # If there are no unclassified points left, stop clustering
    if (!length(unclass))
      break
  }

  
# Create result list and assign class
  result <- list(cluster = cv, eps = eps, 
              eps2 = eps2, eps3 = eps3,
              minpts = minpts, density = classn)
  rm(classn)

  class(result) <- "mdst-dbscan"
  return(result)
}
