#######------------------------------------ Source code for running M-HL codes ------------------------------------#######
comp_HL_gamma1 <- function(formula, EVENT, init, seed){ 
set.seed(seed)

###----- fixed, random effect parameter -----###
HL_fix = function(para,leta) {
	Beta1 = para[1:p]
	Beta2 = para[(p+1):(2*p)]
	g1_vec = exp(para[(2*p+1):(2*p+5)])
	g2_vec = exp(para[(2*p+6):(2*p+10)])
  
	eta = exp(leta)

#--- frailty (ui) ---#
	m_i = function(g) {

		return(sum(event1[group == g]))
	}

	L_i = function(g) {
		bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector( Beta1%*%Z[group == g,] ) } else {
		as.vector( Beta1%*%t(Z[group == g,]) ) }

		return(sum(I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1)))
	}
  
	M_i = function(g) {

		return(sum(event2[group == g]))
	}

	R_i = function(g) {
		bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector( Beta2%*%Z[group == g,] ) } else {
		as.vector( Beta2%*%t(Z[group == g,]) ) }

		return(sum(I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2)))
	}
  
	v_temp1 = sapply(c(1:G),m_i)+sapply(c(1:G),M_i)+1/eta
	v_temp2 = sapply(c(1:G),L_i)+sapply(c(1:G),R_i)+1/eta
	v_vec = log(v_temp1)-log(v_temp2)

	cond_LH_func = function(g) {
		bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta1%*%Z[group == g,]) } else {
		as.vector(Beta1%*%t(Z[group == g,])) }

		l_func = function(t) {M.spline(t,t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])}
		L_func = function(t) {I.spline(t,t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])}
    
		bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta2%*%Z[group == g,]) } else {
		as.vector(Beta2%*%t(Z[group == g,])) }

		r_func = function(t) {M.spline(t,t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])}
		R_func = function(t) {I.spline(t,t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])}
    
		temp1 = event1[group == g]*log(l_func(t_event1[group == g]))-L_func(t_event1[group == g])
		temp2 = event2[group == g]*log(r_func(t_event2[group == g]))-R_func(t_event2[group == g])
		temp = sum(temp1+temp2)
    
		return(temp)
	}

	cond_LH = sapply(c(1:G),cond_LH_func) # ell_1ij
	logfv = (v_vec-exp(v_vec))/eta-lgamma(1/eta)-log(eta)/eta # ell_2ij
  
	h = sum(cond_LH)+sum(logfv)
  
	Omega = c(192,-132,24,12,0,-132,96,-24,-12,12,24,-24,24,-24,24,12,-12,-24,96,-132,0,12,24,-132,192) # penalty term
  
	Omega1 = matrix(Omega,5,5)/((t1max-t1min)/2)^5
	Pen1 = t(g1_vec)%*%Omega1%*%g1_vec
	Omega2 = matrix(Omega,5,5)/((t2max-t2min)/2)^5
	Pen2 = t(g2_vec)%*%Omega2%*%g2_vec

	return(h-kappa1*Pen1-kappa2*Pen2)
} # end of HL_fix

###----- dispersin parameter ----- ###
HL_dis = function(para,leta) {
	Beta1 = para[1:p]
	Beta2 = para[(p+1):(2*p)]
	g1_vec = exp(para[(2*p+1):(2*p+5)])
	g2_vec = exp(para[(2*p+6):(2*p+10)])
  
	eta = exp(leta)

	m_i = function(g) {

		return(sum(event1[group == g]))
	}

	L_i = function(g) {
		bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta1%*%Z[group == g,]) } else {
		as.vector(Beta1%*%t(Z[group == g,])) }

		return(sum(I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1)))
	}
  
	M_i = function(g) {
	
		return(sum(event2[group == g]))
	}

	R_i = function(g) {
		bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta2%*%Z[group == g,]) } else {
		as.vector(Beta2%*%t(Z[group == g,])) }

		return(sum(I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2)))
	}
  
	v_temp1 = sapply(c(1:G),m_i)+sapply(c(1:G),M_i)+1/eta
	v_temp2 = sapply(c(1:G),L_i)+sapply(c(1:G),R_i)+1/eta
	v_vec = log(v_temp1)-log(v_temp2)
  
	cond_LH_func = function(g) {
		bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta1%*%Z[group == g,]) } else {
		as.vector(Beta1%*%t(Z[group == g,])) }

		l_func = function(t) {M.spline(t,t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])}
		L_func = function(t) {I.spline(t,t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])}
    
		bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta2%*%Z[group == g,]) } else {
		as.vector(Beta2%*%t(Z[group == g,])) }

		r_func = function(t) {M.spline(t,t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])}
		R_func = function(t) {I.spline(t,t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])}
    
		temp1 = event1[group == g]*log(l_func(t_event1[group == g]))-L_func(t_event1[group == g])
		temp2 = event2[group == g]*log(r_func(t_event2[group == g]))-R_func(t_event2[group == g])
		temp = sum(temp1+temp2)
    
		return(temp)
	}

	cond_LH = sapply(c(1:G),cond_LH_func) # ell_1ij
	logfv = (v_vec-exp(v_vec))/eta-lgamma(1/eta)-log(eta)/eta # ell_2ij
  
	h = sum(cond_LH)+sum(logfv)

	negative_hess_func = function(g) {
		bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta1%*%Z[group == g,]) } else {
		as.vector(Beta1%*%t(Z[group == g,])) }
		bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(Beta2%*%Z[group == g,]) } else {
		as.vector(Beta2%*%t(Z[group == g,])) }

		comp1 = sum(I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g]))
		comp2 = sum(I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g]))
		temp = comp1+comp2+exp(v_vec[g])/eta
    
		return(temp)
	}

	negative_hess = sapply(c(1:G),negative_hess_func)
  
	adj1 = sum(log(negative_hess/(2*pi))) # Second term of pv(h) in P8 II) Estimation of eta
	adj2 = -2*sum((exp(v_vec)*(sapply(c(1:G),L_i)+sapply(c(1:G),R_i)+1/eta))^-1) # F(h) in P9

	return(h-adj1/2-adj2/24)
} # end of HL_dis

### Simulation setting
eps = 1e-6

Iter_max = 200
Rand_max = 200

p <- dim(X)[2]-2

t_event1 <- EVENT[,1]
t_event2 <- EVENT[,2]
event1 <- EVENT[,3]
event2 <- EVENT[,4]

t1max <- max(t_event1)
t1min <- min(t_event1)
t2max <- max(t_event2)
t2min <- min(t_event2)
  
para_old <- init[1:(length(init)-3)]
leta_ini <- init[length(init)-2]
leta_old <- leta_ini
kappa1 <- init[length(init)-1]
kappa2 <- init[length(init)]

rand = 0

### Estimating results
repeat {
	if (rand >= Rand_max) {break}

	### fixed, random
	res_fix = optim(para_old, HL_fix, leta = leta_old, method = "BFGS",
			control = list(fnscale = -1), hessian = TRUE)
	para_new = res_fix$par

	# singular test
	fix_singular = try(solve(-res_fix$hessian),silent = TRUE)
    
	if (is(fix_singular,"try-error")) {
		leta_old = leta_ini+runif(1,0,2)
		next
	}

	### dispersion
	res_dis = try(optim(leta_old, HL_dis, para = para_new, method = "BFGS",
                        control = list(fnscale = -1),hessian = TRUE),silent = TRUE)
	leta_new = res_dis$par
	sv_hPL = res_dis$value

	# singular test
	if (is(res_dis,"try-error") | solve(-res_dis$hessian) < 0) {
		leta_old = leta_ini+runif(1,0,2)
		next
	}

	dif_old = c(para_old[1:(2*p)],exp(para_old[(2*p+1):(2*p+10)]),exp(leta_old))
	dif_new = c(para_new[1:(2*p)],exp(para_new[(2*p+1):(2*p+10)]),exp(leta_new))
	dif = dif_new-dif_old

	if (max(abs(dif)) < eps) {break}

	para_old = para_new
	leta_old = leta_new

      rand = rand+1
	print(rand)
}

### Results
res_Beta1 <- c(
	para_new[1:p]
	)

res_Beta2 <- c(
	para_new[(p+1):(2*p)]
	)

res_g1 <- rbind(
	exp( para_new[(2*p+1):(2*p+5)] )
	)

res_g2 <- rbind(
	exp( para_new[(2*p+6):(2*p+10)] )
	)

res_eta <- cbind(
	exp(leta_new)
	)


### Computing SE
V <- solve(-res_fix$hessian)
Beta1_se <- sqrt(diag(V)[1:p])
Beta2_se <- sqrt(diag(V)[(p+1):(2*p)])
leta_se <- sqrt(solve(-res_dis$hessian))
eta_se <- res_eta*leta_se

### random effect
m_i = function(g) {
	return(sum(event1[group == g]))
}

L_i = function(g) {
	bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
	as.vector( res_Beta1%*%Z[group == g,] ) } else {
	as.vector( res_Beta1%*%t(Z[group == g,]) ) }

	return( sum( I.spline(t_event1[group == g],t1min,t1max)%*%t(res_g1)*exp(bZ1) ) )
}
  
M_i = function(g) {
	return(sum(event2[group == g]))
}

R_i = function(g) {
	bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
	as.vector( res_Beta2%*%Z[group == g,] ) } else {
	as.vector( res_Beta2%*%t(Z[group == g,]) ) }

	return(sum(I.spline(t_event2[group == g],t2min,t2max)%*%t(res_g2)*exp(bZ2)))
}
  
v_temp1 = sapply(c(1:G),m_i)+sapply(c(1:G),M_i)+1/res_eta
v_temp2 = sapply(c(1:G),L_i)+sapply(c(1:G),R_i)+1/res_eta
v_vec_est = log(v_temp1)-log(v_temp2)
u_vec_est = exp(v_vec_est)

hess_func = function(para, leta) {
	Beta1 = para[1:p]
	Beta2 = para[(p+1):(2*p)]
	g1_vec = exp( para[(2*p+1):(2*p+5)] )
	g2_vec = exp( para[(2*p+6):(2*p+10)] )
	v_vec = para[(2*p+11):(2*p+G+10)]
    
	eta = exp(leta)
    
	m_i = function(g) {
		return(sum(event1[group == g]))
	}

	L_i = function(g) {
		bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
		return(sum(I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1)))
	}
    
	M_i = function(g) {
		return(sum(event2[group == g]))
	}

	R_i = function(g) {
		bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
		return(sum(I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2)))
	}
    
	Mat = matrix(0,2*p+G+10,2*p+G+10)
	Mat_s = matrix(0,2*p+G+10,2*p+G+10)
    
	### Beta1 HESS
	for (i in 1:p) {
		Beta1_diag_func = function(g) {
			bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
			temp = -Z[group == g,i]^2*I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])
        
			return( sum(temp) )
      	}
      Mat[i,i] = Mat_s[i,i] = sum( sapply(c(1:G), Beta1_diag_func) )
	}

	if (p > 1) {
		for (i in 2:p) {
			for (j in 1:(i-1)) {
				Beta1_offdiag_func = function(g) {
					bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
					temp = -Z[group == g,i]*Z[group == g,j]*I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])
            
					return( sum(temp) )
				}

				Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = sum(sapply(c(1:G),Beta1_offdiag_func))
			}
		}
	}
    
	### Beta2 HESS
	for (i in (p+1):(2*p)) {
		Beta2_diag_func = function(g) {
			bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
			temp = -Z[group == g,i-p]^2*I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])
	
			return( sum(temp) )
      	}

	Mat[i,i] = Mat_s[i,i] = sum(sapply(c(1:G),Beta2_diag_func))
	}

	if (p > 1) {
		for (i in (p+2):(2*p)) {
			for (j in (p+1):(i-1)) {
				Beta2_offdiag_func = function(g) {
					bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
					temp = -Z[group == g,i-p]*Z[group == g,j-p]*I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])
            
					return( sum(temp) )
				}

				Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = sum(sapply(c(1:G),Beta2_offdiag_func))
			}
		}
	}
    
	### g1_vec HESS - spline parameter
	for (i in (2*p+1):(2*p+5)) {
		g1_vec1_func = function(g) {
			M.spline_vec = M.spline(t_event1[group == g],t1min,t1max)
			I.spline_vec = I.spline(t_event1[group == g],t1min,t1max)
			bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
        
			comp1 = g1_vec[i-2*p]*M.spline_vec[,i-2*p]/M.spline_vec%*%g1_vec
			comp2 = (g1_vec[i-2*p]*M.spline_vec[,i-2*p]/M.spline_vec%*%g1_vec)^2
			temp = event1[group == g]*(comp1-comp2)-g1_vec[i-2*p]*I.spline_vec[,i-2*p]*exp(bZ1+v_vec[g])
        
			return( sum(temp) )
		}

		Mat[i,i] = Mat_s[i,i] = sum(sapply(c(1:G),g1_vec1_func))
	}

	for (i in (2*p+2):(2*p+5)) {
		for (j in (2*p+1):(i-1)) {
			g1_vec2_func = function(g) {
				M.spline_vec = M.spline(t_event1[group == g],t1min,t1max)
				temp = -event1[group == g]*g1_vec[i-2*p]*M.spline_vec[,i-2*p]*g1_vec[j-2*p]*M.spline_vec[,j-2*p]/(M.spline_vec%*%g1_vec)^2
	
				return( sum(temp) )
			}
	
		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = sum(sapply(c(1:G),g1_vec2_func))
		}
	}
    
	### g2_vec HESS - spline parameter
	for (i in (2*p+6):(2*p+10)) {
		g2_vec1_func = function(g) {
			M.spline_vec = M.spline(t_event2[group == g],t2min,t2max)
			I.spline_vec = I.spline(t_event2[group == g],t2min,t2max)
			bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
        
			comp1 = g2_vec[i-2*p-5]*M.spline_vec[,i-2*p-5]/M.spline_vec%*%g2_vec
			comp2 = (g2_vec[i-2*p-5]*M.spline_vec[,i-2*p-5]/M.spline_vec%*%g2_vec)^2
			temp = event2[group == g]*(comp1-comp2)-g2_vec[i-2*p-5]*I.spline_vec[,i-2*p-5]*exp(bZ2+v_vec[g])
        
			return( sum(temp) )
   	   }

	Mat[i,i] = Mat_s[i,i] = sum(sapply(c(1:G),g2_vec1_func))
	}

	for (i in (2*p+7):(2*p+10)) {
		for (j in (2*p+6):(i-1)) {
			g2_vec2_func = function(g) {
				M.spline_vec = M.spline(t_event2[group == g],t2min,t2max)
				temp = -event2[group == g]*g2_vec[i-2*p-5]*M.spline_vec[,i-2*p-5]*g2_vec[j-2*p-5]*M.spline_vec[,j-2*p-5]/(M.spline_vec%*%g2_vec)^2
          
				return( sum(temp) )
			}

		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = sum(sapply(c(1:G),g2_vec2_func))
		}
	}
    
	### v_vec HESS - random effect
	for (i in (2*p+11):(2*p+G+10)) {
		v_vec1_func = function(g) {
			bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
			bZ2 = as.vector(Beta2%*%t(Z[group == g,]))

			comp1 = -sum(I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g]))
			comp2 = -sum(I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g]))
			temp = comp1+comp2-exp(v_vec[g])/eta
        
			return( sum(temp) )
		}

	Mat[i,i] = v_vec1_func(i-2*p-10)
	}

	### v_vec HESS - random effect - for H* in trace(H-1 H*) in cAIC
	for (i in (2*p+11):(2*p+G+10)) {
		v_vec1_func = function(g) {
			bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
			bZ2 = as.vector(Beta2%*%t(Z[group == g,]))

			comp1 = -sum(I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g]))
			comp2 = -sum(I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g]))
			temp = comp1+comp2+0.01 # to avoid singular issue
        
			return( sum(temp) )
		}

	Mat_s[i,i] = v_vec1_func(i-2*p-10)
	}

	### Beta1 & g1_vec HESS
	for (i in 1:p) {
		for (j in (2*p+1):(2*p+5)) {
			Beta1_g1_vec_func = function(g) {
				I.spline_vec = I.spline(t_event1[group == g],t1min,t1max)
				bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
				temp = -Z[group == g,i]*g1_vec[j-2*p]*I.spline_vec[,j-2*p]*exp(bZ1+v_vec[g])
         
				return( sum(temp) )
			}
			Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = sum(sapply(c(1:G),Beta1_g1_vec_func))
		}
	}
    
	### Beta2 & g2_vec HESS
	for (i in (p+1):(2*p)) {
		for (j in (2*p+6):(2*p+10)) {
			Beta2_g2_vec_func = function(g) {
				I.spline_vec = I.spline(t_event2[group == g],t2min,t2max)
				bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
				temp = -Z[group == g,i-p]*g2_vec[j-2*p-5]*I.spline_vec[,j-2*p-5]*exp(bZ2+v_vec[g])

				return( sum(temp) )
			}
	
		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = sum(sapply(c(1:G),Beta2_g2_vec_func))
		}
	}
    
	### Beta1 & v_vec HESS
	for (i in 1:p) {
		for (j in (2*p+11):(2*p+G+10)) {
			Beta1_vec_func = function(g) {
				bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
				temp = -sum(Z[group == g,i]*I.spline(t_event1[group == g],t1min,t1max)%*%g1_vec*exp(bZ1+v_vec[g])) 
          
				return( sum(temp) )
			}

		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = Beta1_vec_func(j-2*p-10)
		}
	}
    
	### Beta2 & v_vec HESS
	for (i in (p+1):(2*p)) {
		for (j in (2*p+11):(2*p+G+10)) {
			Beta2_vec_func = function(g) {
				bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
				temp = -sum(Z[group == g,i-p]*I.spline(t_event2[group == g],t2min,t2max)%*%g2_vec*exp(bZ2+v_vec[g])) 
          
				return( sum(temp) )
			}

		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = Beta2_vec_func(j-2*p-10)
		}
	}
    
	### g1_vec & v_vec HESS
	for (i in (2*p+1):(2*p+5)) {
  	    for (j in (2*p+11):(2*p+G+10)) {
			g1_vec_v_vec_func = function(g) {
				I.spline_vec = I.spline(t_event1[group == g],t1min,t1max)
				bZ1 = as.vector(Beta1%*%t(Z[group == g,]))
				temp = -sum(g1_vec[i-2*p]*I.spline_vec[,i-2*p]*exp(bZ1+v_vec[g]))

				return( sum(temp) )
			}

		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = g1_vec_v_vec_func(j-2*p-10)
		}
	}
    
	### g2_vec & v_vec HESS
	for (i in (2*p+6):(2*p+10)) {
    	  for (j in (2*p+11):(2*p+G+10)) {
			g2_vec_v_vec_func = function(g) {
				I.spline_vec = I.spline(t_event2[group == g],t2min,t2max)
				bZ2 = as.vector(Beta2%*%t(Z[group == g,]))
				temp = -sum(g2_vec[i-2*p-5]*I.spline_vec[,i-2*p-5]*exp(bZ2+v_vec[g]))

				return( sum(temp) )
			}

		Mat[i,j] = Mat[j,i] = Mat_s[i,j] = Mat_s[j,i] = g2_vec_v_vec_func(j-2*p-10)
      	}
	}
    
	f = function(para) {
		g1_vec = exp(para[1:5])
		g2_vec = exp(para[6:10])
      
      	Omega = c(192,-132,24,12,0,-132,96,-24,-12,12,24,-24,24,-24,24,12,-12,-24,96,-132,0,12,24,-132,192)
      
      	Omega1 = matrix(Omega,5,5)/((t1max-t1min)/2)^5
      	Pen1 = t(g1_vec)%*%Omega1%*%g1_vec
      	Omega2 = matrix(Omega,5,5)/((t2max-t2min)/2)^5
      	Pen2 = t(g2_vec)%*%Omega2%*%g2_vec
      
      	Pen = kappa1*Pen1+kappa2*Pen2
      
      	return(Pen)
	}

	Pen_hess = hessian(f,para[(2*p+1):(2*p+10)])
    
	Mat[(2*p+1):(2*p+10),(2*p+1):(2*p+10)] = Mat[(2*p+1):(2*p+10),(2*p+1):(2*p+10)]-Pen_hess
	Mat_s[(2*p+1):(2*p+10),(2*p+1):(2*p+10)] = Mat_s[(2*p+1):(2*p+10),(2*p+1):(2*p+10)]-Pen_hess

	result = list(a = Mat, b = Mat_s)
	#result = list(Mat)

	return(result)
} # end of hess_func

res_list = hess_func(c(res_Beta1, res_Beta2, log(res_g1), log(res_g2), v_vec_est), log(res_eta) )
V_full = solve( -res_list$a )
V_full_s = -res_list$b
v_vec_se = sqrt(diag(V_full)[(2*p+11):(2*p+G+10)])
u_vec_se = u_vec_est*v_vec_se

### compute ell_1ij for cAIC ###
	cond_LH_func = function(g) {
		bZ1 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(res_Beta1%*%Z[group == g,]) } else {
		as.vector(res_Beta1%*%t(Z[group == g,])) }

		l_func = function(t) {M.spline(t,t1min,t1max)%*%t(res_g1)*exp(bZ1+v_vec_est[g])}
		L_func = function(t) {I.spline(t,t1min,t1max)%*%t(res_g1)*exp(bZ1+v_vec_est[g])}
    
		bZ2 = if ( is.null( dim(Z[group == g,]) ) ) {
		as.vector(res_Beta2%*%Z[group == g,]) } else {
		as.vector(res_Beta2%*%t(Z[group == g,])) }

		r_func = function(t) {M.spline(t,t2min,t2max)%*%t(res_g2)*exp(bZ2+v_vec_est[g])}
		R_func = function(t) {I.spline(t,t2min,t2max)%*%t(res_g2)*exp(bZ2+v_vec_est[g])}
    
		temp1 = event1[group == g]*log(l_func(t_event1[group == g]))-L_func(t_event1[group == g])
		temp2 = event2[group == g]*log(r_func(t_event2[group == g]))-R_func(t_event2[group == g])
		temp = sum(temp1+temp2) # For list (bladder data)
    
		return(temp)
	}

	ell_1ij = sum( sapply(c(1:G), cond_LH_func) ) # ell_1ij
	cAIC = -2*ell_1ij+2*sum( diag( V_full%*%V_full_s ) )

result <- list(round(cbind(res_Beta1, Beta1_se), 3), round(cbind(res_Beta2, Beta2_se), 3),
			round( res_g1, 3 ), round( res_g2, 3 ), round(cbind(res_eta, eta_se), 3)
			, round(sv_hPL, 3), round(-2*sv_hPL+(p+p+5+5+1), 3), round(cAIC, 3), v_vec_est, v_vec_se)

names(result)[[1]]<-"Est_beta 1"
names(result)[[2]]<-"Est_beta 2"
names(result)[[3]]<-"Est_g1"
names(result)[[4]]<-"Est_g2"
names(result)[[5]]<-"Est_eta"
names(result)[[6]]<-"Sv(h)"
names(result)[[7]]<-"mAIC"
names(result)[[8]]<-"cAIC"
names(result)[[9]]<-"Est Random effect"
names(result)[[10]]<-"SE of Random effect"

colnames(result[[1]]) <- c("Est","SE")
colnames(result[[2]]) <- c("Est","SE")
colnames(result[[5]]) <- c("Est","SE")

return(result)

} # end of source



