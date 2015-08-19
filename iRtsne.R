iRtsne <- function(X, pca=TRUE, dims=2, initial_dims=NULL, perplexity=30, walkmap = FALSE,
                   seed=NULL, max_iter=1000, epoch=100, min_cost=0.5 ){
    ## initialize variables
    if (is.null(initial_dims) | missing(initial_dims)) 
        initial_dims <- ncol(X)
    if(!(is.null(seed))) set.seed(seed)          # set a seed to fix the initial embedding
    n = nrow(X)
    momentum = 0.5
    final_momentum = 0.8
    mom_switch_iter = 250
    epsilon = 500
    min_gain = 0.01
    initial_P_gain = 4                           # early exaggeration of P
    eps = 2^(-52)                                # typical machine precision
    ydata = matrix(rnorm(dims*n), n, dims)       # initialize map results data
    grads =  matrix(0,nrow(ydata),ncol(ydata))   # gradient descent
    incs =  matrix(0,nrow(ydata),ncol(ydata))
    gains = matrix(1,nrow(ydata),ncol(ydata))
    
    ## Apply PCA to tidy the data
    X = as.matrix(X)
    if (pca) {
        cat('Apply PCA to tidy the data and initialize the embedding...')
        pca_result <- prcomp(X, retx=TRUE)
        X <- pca_result$x[ ,1:min(initial_dims,ncol(pca_result$x))]
        #ydata <- X[ ,1:dims]
        cat("DONE!\n")
    }
    
	## P calculation
    if(walkmap == FALSE){
        x2p = calx2p(X, perplexity, 1e-5)
        P <- x2p$P
        sigma <- x2p$sigma
        P = (P + t(P)) / (2*sum(P))
    }else{
        P <- walkmap(X, pca = FALSE, perplexity = perplexity)
        sigma <- 0
    }
	
	P[P < eps] <- eps
	P = P * initial_P_gain
	
	## embedding by minimizing the KL divergence between Q and P
	cat("Optimize embedding by minimizing the KL divergence between Q and P:\n")
	cost <- Inf
	
	library(animation)
	cluster <- c(rep(1, 200), rep(2, 500), rep(3, 1000))
	saveGIF({
	    for (iter in 1:max_iter){
	        if (iter %% epoch == 0) {
	            oldcost <- cost
	            cost <-  sum( P*log((P+eps)/(Q+eps)) )
	            cat(paste0("  Epoch: Iteration #",iter,", KL divergence cost is: ",cost,"\n"))
	            if (cost < min_cost || (oldcost-cost) < 0.01) break  # for speed up
	        }
	        
	        plot(ydata, col=factor(cluster))
	        ani.pause(1)
	        
	        ## calculate Q value  (equal implementation: num<-1/(1+as.matrix(dist(ydata))^2) )
	        sum_ydata = apply(ydata^2, 1, sum)
	        num =  1/(1 + sum_ydata + sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata))) 
	        diag(num)=0
	        Q = num / sum(num)
	        Q[Q < eps] = eps
	        
	        ## Perform gradient update (with momentum and gains)
	        if (any(is.nan(num))) message ('NaN in grad. descent')
	        stiffnesses = 4 * (P-Q) * num
	        for (i in 1:n){
	            grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum) }
	        
	        gains = (gains+0.2)*abs(sign(grads) != sign(incs)) + gains*0.8*abs(sign(grads) == sign(incs))		
	        gains[gains < min_gain] = min_gain
	        incs = momentum * incs - epsilon * (gains * grads)
	        ydata = ydata + incs
	        ydata = sweep(ydata, 2, apply(ydata,2,mean))
	        if (iter == mom_switch_iter) momentum = final_momentum
	        if (iter == 100) P = P/4
	    }
	},  movie.name = "test.gif", ani.width = 800, ani.height = 800)
	
	
# 	for (iter in 1:max_iter){
# 	    if (iter %% epoch == 0) {
# 	        oldcost <- cost
# 	        cost <-  sum( P*log((P+eps)/(Q+eps)) )
# 	        cat(paste0("  Epoch: Iteration #",iter,", KL divergence cost is: ",cost,"\n"))
# 	        if (cost < min_cost || (oldcost-cost) < 0.01) break  # for speed up
# 	    }
# 	    
# 	    ## calculate Q value  (equal implementation: num<-1/(1+as.matrix(dist(ydata))^2) )
# 	    sum_ydata = apply(ydata^2, 1, sum)
# 	    num =  1/(1 + sum_ydata + sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata))) 
# 	    diag(num)=0
# 	    Q = num / sum(num)
# 	    Q[Q < eps] = eps
# 	    
# 	    ## Perform gradient update (with momentum and gains)
# 	    if (any(is.nan(num))) message ('NaN in grad. descent')
# 	    stiffnesses = 4 * (P-Q) * num
# 	    for (i in 1:n){
# 	        grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum) }
# 	    
# 	    gains = (gains+0.2)*abs(sign(grads) != sign(incs)) + gains*0.8*abs(sign(grads) == sign(incs))		
# 	    gains[gains < min_gain] = min_gain
# 	    incs = momentum * incs - epsilon * (gains * grads)
# 	    ydata = ydata + grads  #incs
# 	    ydata = sweep(ydata, 2, apply(ydata,2,mean))
# 	    if (iter == mom_switch_iter) momentum = final_momentum
# 	    if (iter == 100) P = P/4
# 	}
	
	cat("DONE!\n")
	res <- list(Y = ydata, sigma = sigma, P=P, Q=Q)
	class(res) <- 'iRtsne'
	return(res)
}



calx2p <- function(X, perplexity = 30, tol = 1e-4){
    cat("Calculate the pairwise P value...")
	n <- nrow(X)
	D <- as.matrix(dist(X))
	P <- matrix(0, n, n )		
	beta <- rep(1, n)           # beta = 0.5/sigma^2
	## test with variate perplexity
	if(length(perplexity) < n)
	    perplexity <- rep(perplexity, n)
	logU <- log2(perplexity)    # logU equal to entropy H
	
	for (i in 1:n){
		Di <- D[i, -i]
		beta[i] <- ifelse(i>1, beta[i-1], beta[i])
		betamin <- -Inf
		betamax <- Inf
		tries <- 0
		while(tries < 1000){
		    hbeta <- calHbeta(Di, beta[i])
		    H <- hbeta$H
		    thisP <- hbeta$P
		    Hdiff <- H - logU[i]
		    
		    if (abs(Hdiff) <= tol) break
			if (Hdiff > 0){
				betamin <- beta[i]
				if (is.infinite(betamax)){
				    beta[i] <- beta[i] * 2
				}else {
				    beta[i] <- (beta[i] + betamax)/2
				    }
			} else{
				betamax <- beta[i]
				if (is.infinite(betamin)){
				    beta[i] <- beta[i]/ 2
				}else {
				    beta[i] <- (beta[i]+betamin)/2
				    }
				}
			tries <- tries + 1
			}
		P[i,-i] <- thisP	
	}	
	cat("DONE!\n")
	r <- list()
	r$P <- P
	r$sigma <- sqrt(0.5/beta)
	r 
}



## calculate entropy H and probability P with beta
calHbeta <- function(D, beta){
    P = exp(-D^2 * beta)
    eps = 2^(-52) 
    sumP = sum(P)
    if (sumP == 0){
        P = D * 0
        H = 0
    } else {
        P = P/sumP
        P[P < eps] <- eps
        H = -sum(P*log2(P))
    }
    r = list()
    r$P = P
    r$H = H
    r
}















##-----------------------walkmap to reweight P-------------------------------##

walkmap <- function(X, pca = TRUE, perplexity=30){
    
    X = as.matrix(X)
    n <- nrow(X)
    
    ## Apply PCA to tidy the data
    if (pca) {
        cat('Apply PCA to tidy the data and initialize the embedding...')
        pca_result <- prcomp(X, retx=TRUE)
        X <- pca_result$x
    }
    
    ## P calculation, diagnol random walk probability matrix
    x2p = calx2p(X, perplexity, 1e-5)
    P <- x2p$P
    P = (P + t(P)) / 2
    for(i in 1:nrow(P)){
        for (j in 1:ncol(P)){
            if(i>=j) P[i,j] <- 0
        }
    }
    
    ## reweight probability with random walk
    newP <- matrix(0, n, n)    # initialize new P matrix
    thisWalk <- rep(0,n)       # initialize first walk vector
    pb <- txtProgressBar(min = 0, max = n, style = 3) # create progress bar
    for(i in 1:n){
        setTxtProgressBar(pb, i)
        thisWalk[i] <- 1
        while(sum(thisWalk) > 0){
            nextWalk <- randomWalk(thisWalk, P)
            thisWalk <- nextWalk
            newP[,i] <- newP[,i] + nextWalk
        }
    }
    close(pb)
    
    newP <- t(newP)
    for(i in 1:nrow(newP)){
        for (j in 1:ncol(newP)){
            if(i>=j) newP[i,j] <- newP[j,i]
        }
    }
    
    P = newP/sum(newP) 
    return(P)
}


randomWalk <- function(v, m){
    o <- v %*% m
    return(t(o)[,1])
}


## conditional probability between point i and j.
## for Pj|i, neighbors are picked in proportion to their probability density(perplexity) 
## under a gaussian centered at point i 

calx2p <- function(X, perplexity = 30, tol = 1e-4){
    cat("Calculate the pairwise P value...")
    n <- nrow(X)
    D <- as.matrix(dist(X))
    P <- matrix(0, n, n )		
    beta <- rep(1, n)           # beta = 0.5/sigma^2
    ## test with variate perplexity
    if(length(perplexity) < n)
        perplexity <- rep(perplexity, n)
    logU <- log2(perplexity)    # logU equal to entropy H
    
    for (i in 1:n){
        Di <- D[i, -i]
        beta[i] <- ifelse(i>1, beta[i-1], beta[i])
        betamin <- -Inf
        betamax <- Inf
        tries <- 0
        while(tries < 1000){
            hbeta <- calHbeta(Di, beta[i])
            H <- hbeta$H
            thisP <- hbeta$P
            Hdiff <- H - logU[i]
            
            if (abs(Hdiff) <= tol) break
            if (Hdiff > 0){
                betamin <- beta[i]
                if (is.infinite(betamax)){
                    beta[i] <- beta[i] * 2
                }else {
                    beta[i] <- (beta[i] + betamax)/2
                }
            } else{
                betamax <- beta[i]
                if (is.infinite(betamin)){
                    beta[i] <- beta[i]/ 2
                }else {
                    beta[i] <- (beta[i]+betamin)/2
                }
            }
            tries <- tries + 1
        }
        P[i,-i] <- thisP	
    }	
    cat("DONE!\n")
    r <- list()
    r$P <- P
    r$sigma <- sqrt(0.5/beta)
    r 
}



## calculate entropy H and probability P with beta
calHbeta <- function(D, beta){
    P = exp(-D^2 * beta)
    eps = 2^(-52) 
    sumP = sum(P)
    if (sumP == 0){
        P = D * 0
        H = 0
    } else {
        P = P/sumP
        P[P < eps] <- eps
        H = -sum(P*log2(P))
    }
    P[P <= eps] <- 0
    r = list()
    r$P = P
    r$H = H
    r
}
