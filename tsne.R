tsne <- function(X, pca = TRUE, dims=2, initial_dims=NULL, perplexity=30, 
                 max_iter = 1000, min_cost=0, epoch_callback=NULL,whiten=TRUE, epoch=100 ){
    
    ## initialize variables
    momentum = 0.5
    final_momentum = 0.8
    mom_switch_iter = 250
    epsilon = 500
    min_gain = 0.01
    initial_P_gain = 4
    eps = 2^(-52)                  # typical machine precision
    
    ## Apply PCA
    X = as.matrix(X)
    if (pca) {
        pca_result <- prcomp(X,retx=TRUE)
        X <- pca_result$x[,1:min(initial_dims,ncol(pca_result$x))]
    }
    
    ## initialize map data
    if (is.null(initial_dims) | missing(initial_dims)) 
        initial_dims <- ncol(X)
    n = nrow(X)
	ydata = matrix(rnorm(dims * n),n)
	
	## P calculation
	P = .x2p(X, perplexity, 1e-5)$P
	P = .5 * (P + t(P))
	P[P < eps]<-eps
	P = P/sum(P)
	
	
	P = P * initial_P_gain
	grads =  matrix(0,nrow(ydata),ncol(ydata))
	incs =  matrix(0,nrow(ydata),ncol(ydata))
	gains = matrix(1,nrow(ydata),ncol(ydata))

	
	for (iter in 1:max_iter){
		if (iter %% epoch == 0) { # epoch
			cost =  sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
			message("Epoch: Iteration #",iter," error is: ",cost)
			if (cost < min_cost) break
			if (!is.null(epoch_callback)) epoch_callback(ydata)
		}


		sum_ydata = apply(ydata^2, 1, sum)
		num =  1/(1 + sum_ydata +    sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata))) 
		diag(num)=0
		Q = num / sum(num)
		if (any(is.nan(num))) message ('NaN in grad. descent')
		Q[Q < eps] = eps
		stiffnesses = 4 * (P-Q) * num
		for (i in 1:n){
			grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum)
		}
		
		gains = (gains + .2) * abs(sign(grads) != sign(incs)) 
				+ gains * .8 * abs(sign(grads) == sign(incs))		
		gains[gains < min_gain] = min_gain

		incs = momentum * incs - epsilon * (gains * grads)
		ydata = ydata + incs
		ydata = sweep(ydata,2,apply(ydata,2,mean))
		if (iter == mom_switch_iter) momentum = final_momentum
		
		if (iter == 100) P = P/4
	}
	ydata
}




.Hbeta <- function(D, beta){
	P = exp(-D^2 * beta)
	sumP = sum(P)
	if (sumP == 0){
	    P = D * 0
		H = 0
	} else {
	    P = P/sumP
		H = - sum(P*log2(P))
	}
	r = list()
	r$P = P
	r$H = H
	r
}

.x2p <- function(X, perplexity = 30, tol = 1e-5){

	n = nrow(X)
	D = as.matrix(dist(X))
	P = matrix(0, n, n )		
	beta = rep(1, n)           # beta = 0.5/sigma^2
	logU = log2(perplexity)    # logU equal to entropy H
	
	for (i in 1:n){
		betamin = -Inf
		betamax = Inf
		Di = D[i, -i]
		hbeta = .Hbeta(Di, beta[i])
		H = hbeta$H
		thisP = hbeta$P
		Hdiff = H - logU
		tries = 0

		while(abs(Hdiff) > tol && tries < 50){
			if (Hdiff > 0){
				betamin = beta[i]
				if (is.infinite(betamax)) beta[i] = beta[i] * 2
				else beta[i] = (beta[i] + betamax)/2
			} else{
				betamax = beta[i]
				if (is.infinite(betamin))  beta[i] = beta[i]/ 2
				else beta[i] = ( beta[i] + betamin) / 2
			}
			
			hbeta = .Hbeta(Di, beta[i])
			H = hbeta$H
			thisP = hbeta$P
			Hdiff = H - logU
			tries = tries + 1
		}	
			P[i,-i]  = thisP	
	}	
	
	r = {}
	r$P = P
	r$beta = beta
	
	sigma = sqrt(0.5/beta)
	message('sigma summary: ', paste(names(summary(sigma)),':',summary(sigma),'|',collapse=''))

	r 
}

.whiten <- function(X, row.norm=FALSE, verbose=FALSE, n.comp=ncol(X))
{  
	n.comp; # forces an eval/save of n.comp
	if (verbose) message("Centering")
   n = nrow(X)
	p = ncol(X)
	X <- scale(X, scale = FALSE)
   X <- if (row.norm) 
       t(scale(X, scale = row.norm))
   else t(X)

   if (verbose) message("Whitening")
   V <- X %*% t(X)/n
   s <- La.svd(V)
   D <- diag(c(1/sqrt(s$d)))
   K <- D %*% t(s$u)
   K <- matrix(K[1:n.comp, ], n.comp, p)
   X = t(K %*% X)
	X
}

