# This is a modified form of the IMIS function from package "IMIS"
# Parallelized version, and customized to this model

IMIS_mod_par <- function (B = 1000, B.re = 3000, number_k = 100, D = 0) {
	
	## Setting up parallelization
	# # Allows up to 8 runs in parallel - customize to your machine
	max_clusters <- 8
	cl1 <- makeCluster(max_clusters)
	junk <- clusterEvalQ(cl1, library(deSolve))
	junk2 <- clusterEvalQ(cl1, library(logitnorm))
	clusterExport(cl1, c("malaria_model", "calc_monthly_west_llik", "calc_monthly_east_llik", "calc_prev_west_llik", "calc_prev_east_llik", "calc_pfmdr_west_llik", "calc_pfmdr_east_llik", "calc_pfpm2_west_llik", "calc_pfpm2_east_llik", "calc_mtes_west_llik", "calc_mtes_east_llik", "calc_ptes_west_llik", "calc_ptes_east_llik", 
	"prior","llikelihood","rmvnorm")) #just added these 2
	clusterExport(cl1, c("B"),envir=environment())
	
	# set up for initial run
    B0 = B * 10
    X_all = X_k = sample.prior(B0)
    X_all = as.matrix(X_all)
    if (is.vector(X_all)) 
        Sig2_global = var(X_all)
    if (is.matrix(X_all)) 
        Sig2_global = cov(X_all)
    stat_all = matrix(NA, 6, number_k)
    center_all = prior_all = like_all = llike_all = NULL
    sigma_all = list()
    if (D >= 1) 
        option.opt = 1
    if (D == 0) {
        option.opt = 0
        D = 1
    }
    for (k in 1:number_k) { # number of IMIS iterations
    		print('hello1') #to keep track of location
        ptm.like = proc.time()
        X_k <- as.data.frame(X_k)
        prior_tmp <- prior(X_k)
        prior_all = c(prior_all, prior_tmp)
       
        tmp_llik <- parApply(cl1, X_k, 1, llikelihood)
        tmp_llik[prior_tmp==0] <- NA 
        llike_all = c(llike_all, tmp_llik)
        llike_all[is.nan(llike_all)] <- -Inf
        llike_all[is.na(llike_all)] <- -Inf
        print('log likelihood')
        print(sort(llike_all))
         important_ak = which(llike_all == max(llike_all))
 		X_imp_ak = X_all[important_ak, ]
        like_all = exp(llike_all - max(llike_all))
        ptm.use = (proc.time() - ptm.like)[3]
        if (ptm.use<1) { error}
        if (k == 1) 
            print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60, 
                2), "minutes"))
        if (k == 1) 
            envelop_all = prior_all
        if (k > 1) 
            envelop_all = apply(rbind(prior_all * B0/B, gaussian_all), 
                2, sum)/(B0/B + D + (k - 2))               
        Weights = prior_all * like_all/envelop_all
        stat_all[1, k] = log(mean(Weights))
        Weights = Weights/sum(Weights)
        stat_all[2, k] = sum(1 - (1 - Weights)^B.re)
        stat_all[3, k] = max(Weights)
        stat_all[4, k] = 1/sum(Weights^2)
        stat_all[5, k] = -sum(Weights * log(Weights), na.rm = TRUE)/log(length(Weights))
        stat_all[6, k] = var(Weights/mean(Weights))
        if (k == 1) 
            print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
        print(c(k, round(stat_all[1:4, k], 3)))
        print('hello...')
        if (k == 1 & option.opt == 1) {
        	
        	num_opt = D # input how many points to do
        	
             if (is.matrix(X_all)) 
                Sig2_global = cov(X_all[which(llike_all > min(llike_all)), 
                   ])
            X_k = which_exclude = NULL
           label_weight = sort(llike_all, decreasing = TRUE, index = TRUE)
            which_remain = which(llike_all > label_weight$x[num_opt+1]) # select num_opt best points
            																										# as starting points for optimization
            size_remain = length(which_remain)
            print(size_remain)
            X_imp_mat = X_all[which_remain, ] 
			optimizing_imis <- function(X_imp) {
				center_all = prior_all = NULL
				posterior = function(theta) {
                  -log(prior(theta)) - llikelihood(theta)
                }
                
                # Optimize first by Nelder-Mead
				optimizer = optim(X_imp, posterior, method = "Nelder-Mead", 
                    control= c(trace=1, maxit=1000))  
                    print('here')
                  theta.NM = optimizer$par
                  
                    print(theta.NM)
                    
                    #continue optimization with L-BFGS-B
                    # optimization is better when using the 2 methods successively
                	optimizer = optim(theta.NM, posterior, method = "L-BFGS-B",  
                    hessian = TRUE, lower=c(1e-2, 2e-2,1e-2,1e-2,1e-2, 1, 1e-2), 
                    upper=c(1,0.7,1,1,20, 1e6, 1)-1e-2, control = list(trace=6,maxit = 50)) 
                     #ranges for a,  p, f_m, f_p, mosq_m, N_init, beta_min_m
                    
                    print(optimizer$par)
                   
                center_all = optimizer$par
                  if (min(eigen(optimizer$hessian)$values) > 
                    0) 
                    sigma_all_tmp = solve(optimizer$hessian)
                  if (min(eigen(optimizer$hessian)$values) <= 0) {
                    eigen.values = eigen(optimizer$hessian)$values
                    eigen.values[which(eigen.values < 0)] = 0
                    hessian = eigen(optimizer$hessian)$vectors %*% 
                      diag(eigen.values) %*% t(eigen(optimizer$hessian)$vectors)
                    sigma_all_tmp = solve(hessian + diag(1/diag(Sig2_global)))
                  }
                  
                  	X_new = rmvnorm(B, optimizer$par, sigma_all_tmp) 
                  	
                  	return(list(X_new, center_all, sigma_all_tmp))
                  
            }  # end optimizing_imis
                      
            ptm.opt = proc.time()
			
			####### WHAT TO OUTPUT
			print(X_imp_mat)
			 X_new_list <- parApply(cl1, X_imp_mat, 1, optimizing_imis) # run optimizations in parallel as well
        	    X_k_list <- lapply(X_new_list, "[[", 1) 
        	    center_all_list <- lapply(X_new_list, "[[", 2) 
        	    sigma_all <- lapply(X_new_list, "[[", 3) 
        	    
            X_k <- do.call("rbind", X_k_list)
            center_all <-do.call("rbind", center_all_list)
            ptm.use = (proc.time() - ptm.opt)[3]
            print(paste("time used=", round(ptm.use/60, 2)))
                  
            X_all = rbind(X_all, X_k)
                    
        } # end k==1
        if (k > 1 | option.opt == 0) { 
            important = which(Weights == max(Weights))
            if (length(important) > 1) 
                important = important[1]
            if (is.matrix(X_all)) 
                X_imp = X_all[important, ]
            if (is.vector(X_all)) 
                X_imp = X_all[important]
            if (is.matrix(X_all)) 
                center_all = rbind(center_all, X_imp)
            if (is.vector(X_all)) 
                center_all = c(center_all, X_imp)
            if (is.matrix(X_all)) 
                distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)))
            if (is.vector(X_all)) 
                distance_all = abs(X_all - X_imp)
            label_nr = sort(distance_all, decreasing = FALSE, 
                index = TRUE)
            which_var = label_nr$ix[1:B]
            if (is.matrix(X_all)) 
                Sig2 = cov.wt(X_all[which_var, ], wt = Weights[which_var] + 
                  1/length(Weights), cor = FALSE, center = X_imp, 
                  method = "unbias")$cov
            if (is.vector(X_all)) {
                Weights_var = Weights[which_var] + 1/length(X_all)
                Weights_var = Weights_var/sum(Weights_var)
                Sig2 = (X_all[which_var] - X_imp)^2 %*% Weights_var
            }
            sigma_all[[D + k - 1]] = Sig2
            if (is.matrix(X_all)) 
                X_k = rmvnorm(B, X_imp, Sig2)
            if (is.vector(X_all)) 
                X_k = rnorm(B, X_imp, sqrt(Sig2))
            if (is.matrix(X_all)) 
                X_all = rbind(X_all, X_k)
            if (is.vector(X_all)) 
                X_all = c(X_all, X_k)
        }
        if (k == 1) {
        print('hello4')
            gaussian_all = matrix(NA, D, B0 + D * B)
            for (i in 1:D) {
                if (is.matrix(X_all)) 
                  gaussian_all[i, ] = dmvnorm(X_all, center_all[i, 
                    ], sigma_all[[i]])
                if (is.vector(X_all)) 
                  gaussian_all[i, ] = dnorm(X_all, center_all[i], 
                    sqrt(sigma_all[[i]]))
            }
            print('ello')
        }
        if (k > 1) {
        print('ehllo')
            if (is.vector(X_all)) 
                gaussian_new = matrix(0, D + k - 1, length(X_all))
            if (is.matrix(X_all)) 
                gaussian_new = matrix(0, D + k - 1, dim(X_all)[1])
            if (is.matrix(X_all)) {
                gaussian_new[1:(D + k - 2), 1:(dim(X_all)[1] - 
                  B)] = gaussian_all
                gaussian_new[D + k - 1, ] = dmvnorm(X_all, X_imp, 
                  sigma_all[[D + k - 1]])
                for (j in 1:(D + k - 2)) gaussian_new[j, (dim(X_all)[1] - 
                  B + 1):dim(X_all)[1]] = dmvnorm(X_k, center_all[j, 
                  ], sigma_all[[j]])
            }
            if (is.vector(X_all)) {
                gaussian_new[1:(D + k - 2), 1:(length(X_all) - 
                  B)] = gaussian_all
                gaussian_new[D + k - 1, ] = dnorm(X_all, X_imp, 
                  sqrt(sigma_all[[D + k - 1]]))
                for (j in 1:(D + k - 2)) gaussian_new[j, (length(X_all) - 
                  B + 1):length(X_all)] = dnorm(X_k, center_all[j], 
                  sqrt(sigma_all[[j]]))
            }
            gaussian_all = gaussian_new
        }
        if (stat_all[2, k] > (1 - exp(-1)) * B.re) 
            break
    }
    print('hello5')
    nonzero = which(Weights > 0)
    which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
    if (is.matrix(X_all)) 
        resample_X = X_all[which_X, ]
    if (is.vector(X_all)) 
        resample_X = X_all[which_X]
    stopCluster(cl1)
    return(list(stat = t(stat_all), resample = resample_X, center = center_all))
}

