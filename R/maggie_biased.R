#' Test of Synchronization



epanech_kernel<- function(u){
  return(ifelse(abs(u)<=1, 3/4*(1-u^2), 0))
}

S0_gap<- function(i,j,b){
  s<- convolve(rep(1,j-i+1), epanech_kernel((i-j-1):(j-i+1)/(n*b)), type="open")[(j-i+2):(2*(j-i)+2)]
  return(s)
}
S1_gap<- function(i,j,b){
  s<- convolve(rep(1,j-i+1), epanech_kernel((i-j-1):(j-i+1)/(n*b))*((i-j-1):(j-i+1)/(n)),
               type="open")[(j-i+2):(2*(j-i)+2)]
  return(s)
}
S2_gap<- function(i,j,b){
  s<- convolve(rep(1,j-i+1), epanech_kernel((i-j-1):(j-i+1)/(n*b))*((i-j-1):(j-i+1)/(n))^2,
               type="open")[(j-i+2):(2*(j-i)+2)]
  return(s)
}

norm<-function(x){return(sqrt(sum(x^2)))}



cp_test_boot<- function(X, b, nboot=500){
  d<- ncol(X); theta_sup<- c(); tau<-c() ; theta.cum<-rep(0,right_lim+1-left_lim)
  index.cp<-c(); n<- nrow(X); grid<- 1:n/n
  Y<- matrix(nrow=n, ncol=d) ; combn_fit<- matrix(nrow=n, ncol=d)
  combn.2_fit<- matrix(nrow=n, ncol=d) ; combn_fit.og<- matrix(nrow=n, ncol=d)
  S0_mat<- array(dim=c(d, right_lim-left_lim+1,2));S1_mat<- array(dim=c(d, right_lim-left_lim+1,2))
  S2_mat<- array(dim=c(d, right_lim-left_lim+1,2))

  theta_val<-matrix(nrow=right_lim-left_lim+1, ncol=d)

  for(r in 1:d){

    dat_scb<- X[,r]

    S0_mat[r, ,1]<- convolve(rep(1,n), epanech_kernel((2*n*b-1):(-n*b)/(n*b)), type="open")[(left_lim):(right_lim)]
    S0_mat[r, ,2]<- convolve(rep(1,n), epanech_kernel((n*b):(1-2*n*b)/(n*b)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]

    S1_mat[r, ,1]<- convolve(rep(1,n), epanech_kernel((2*n*b-1):(-n*b)/(n*b))*((2*n*b-1):(-n*b)/(n)), type="open")[(left_lim):(right_lim)]
    S1_mat[r, ,2]<- convolve(rep(1,n), epanech_kernel((n*b):(1-2*n*b)/(n*b))*((n*b):(1-2*n*b)/(n)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]

    S2_mat[r, ,1]<- convolve(rep(1,n), epanech_kernel((2*n*b-1):(-n*b)/(n*b))*((2*n*b-1):(-n*b)/(n))^2, type="open")[(left_lim):(right_lim)]
    S2_mat[r, ,2]<- convolve(rep(1,n), epanech_kernel((n*b):(1-2*n*b)/(n*b))*((n*b):(1-2*n*b)/(n))^2, type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]

    f_tilde_L <- convolve(dat_scb[1:n], epanech_kernel((2*n*b-1):(-n*b)/(n*b)), type="open")[(left_lim):(right_lim)]
    f_tilde1_L <- convolve(dat_scb[1:n], epanech_kernel((2*n*b-1):(-n*b)/(n*b))*((2*n*b-1):(-n*b)/(n)), type="open")[(left_lim):(right_lim)]

    dat_fit_L <-   (S2_mat[r, (left_lim:right_lim)-left_lim+1,1]*f_tilde_L -
                      S1_mat[r, (left_lim:right_lim)-left_lim+1,1]*f_tilde1_L)/
      (S2_mat[r, (left_lim:right_lim)-left_lim+1,1]*S0_mat[r, (left_lim:right_lim)-left_lim+1,1]-
         S1_mat[r, (left_lim:right_lim)-left_lim+1,1]^2)

    f_tilde_R <- convolve(dat_scb[1:n], epanech_kernel((n*b):(1-2*n*b)/(n*b)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]
    f_tilde1_R <- convolve(dat_scb[1:n], epanech_kernel((n*b):(1-2*n*b)/(n*b))*((n*b):(1-2*n*b)/(n)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]

    dat_fit_R <-   (S2_mat[r, (left_lim:right_lim)-left_lim+1,2]*f_tilde_R -
                      S1_mat[r, (left_lim:right_lim)-left_lim+1,2]*f_tilde1_R)/
      (S2_mat[r, (left_lim:right_lim)-left_lim+1,2]*S0_mat[r, (left_lim:right_lim)-left_lim+1,2]-
         S1_mat[r, (left_lim:right_lim)-left_lim+1,2]^2)


    dat_fit_L_og <- 2*dat_fit_L - dat_fit_L
    dat_fit_R_og <- 2*dat_fit_R - dat_fit_R
    theta_val[,r]<- abs(dat_fit_L_og-dat_fit_R_og)^2

    theta.cum<- theta_val[,r]+theta.cum
    theta_sup[r] <- max(theta_val[,r])
    tau[r]<- grid[which.max(theta_val[,r])+left_lim-1]
    index.cp[r]<- which.max(theta_val[,r])+left_lim-1;ind<- index.cp[r]

    s0.left<- S0_gap(1, ind, b); s1.left<- S1_gap(1, ind, b); s2.left<- S2_gap(1, ind, b)

    f_tilde.left <- convolve(dat_scb[1:ind], epanech_kernel((-ind):(ind)/(n*b)), type="open")[(ind+1):(2*ind)]
    f_tilde_1.left <-  convolve(dat_scb[1:ind], epanech_kernel((-ind):ind/(n*b))*((-ind):ind/(n)),
                                type="open")[(ind+1):(ind+ind)]
    dat_fit_ll.left <-   (s2.left*f_tilde.left - s1.left * f_tilde_1.left)/(s2.left * s0.left- s1.left^2)

    dat_fit_ll.left.og <- 2*dat_fit_ll.left - dat_fit_ll.left

    s0.right<- S0_gap(ind+1, n,b); s1.right<- S1_gap(ind+1, n, b); s2.right<- S2_gap(ind+1, n, b)
    f_tilde.right <- convolve(dat_scb[(ind+1):n], epanech_kernel((ind-n):(n-ind)/(n*b)), type="open")[(n-ind+1):(2*(n-ind))]
    f_tilde_1.right <-  convolve(dat_scb[(ind+1):n], epanech_kernel((ind-n):(n-ind)/(n*b))*((ind-n):(n-ind)/(n)),
                                 type="open")[(n-ind+1):(2*(n-ind))]
    dat_fit_ll.right <-   (s2.right*f_tilde.right - s1.right * f_tilde_1.right)/(s2.right * s0.right- s1.right^2)


    dat_fit_ll.right.og <- 2*dat_fit_ll.right - dat_fit_ll.right

    s0.combn<- S0_gap(1, n, b); s1.combn<- S1_gap(1, n, b); s2.combn<- S2_gap(1, n, b)
    f_tilde.combn <- convolve(dat_scb[1:n], epanech_kernel((-n):(n)/(n*b)), type="open")[(n+1):(2*n)]
    f_tilde_1.combn <-  convolve(dat_scb[1:n], epanech_kernel((-n):n/(n*b))*((-n):n/(n)),
                                 type="open")[(n+1):(n+n)]
    combn_fit[,r] <-   (s2.combn*f_tilde.combn - s1.combn * f_tilde_1.combn)/(s2.combn * s0.combn- s1.combn^2)

    combn_fit.og[,r]<- 2*combn_fit[,r]- combn_fit[,r]

    Y[,r]<- X[,r] - c(dat_fit_ll.left.og, dat_fit_ll.right.og)

  }
  #et1<-Sys.time()
  theta.combn<- max(theta.cum); tau.combn<- grid[which.max(theta.cum)+left_lim-1]
  ind.combn <- which.max(theta.cum)+left_lim-1
  ts<- abs(theta.combn - sum(theta_sup))

  ### for test 2, fit jointly.
  fit.test.2<- matrix(nrow=n, ncol=d)
  for(r in 1:d){
    # theta.norm.combn<- apply(theta_val, 1, norm)
    # tau.max.combn<- grid[which.max(theta.norm.combn)+left_lim-1]
    dat_scb<- X[,r]
    s0.left.combn<- S0_gap(1, ind.combn, b); s1.left.combn<- S1_gap(1, ind.combn, b); s2.left.combn<- S2_gap(1, ind.combn, b)
    f_tilde.left.combn <- convolve(dat_scb[1:ind.combn], epanech_kernel((-ind.combn):(ind.combn)/(n*b)), type="open")[(ind.combn+1):(2*ind.combn)]
    f_tilde_1.left.combn <-  convolve(dat_scb[1:ind.combn], epanech_kernel((-ind.combn):ind.combn/(n*b))*((-ind.combn):ind.combn/(n)),
                                      type="open")[(ind.combn+1):(ind.combn+ind.combn)]
    dat_fit_ll.left.combn <-   (s2.left.combn*f_tilde.left.combn - s1.left.combn * f_tilde_1.left.combn)/(s2.left.combn * s0.left.combn- s1.left.combn^2)

    s0.right.combn<- S0_gap(ind.combn+1, n, b); s1.right.combn<- S1_gap(ind.combn+1, n, b); s2.right.combn<- S2_gap(ind.combn+1, n, b)
    f_tilde.right.combn <- convolve(dat_scb[(ind.combn+1):n], epanech_kernel((ind.combn-n):(n-ind.combn)/(n*b)), type="open")[(n-ind.combn+1):(2*(n-ind.combn))]
    f_tilde_1.right.combn <-  convolve(dat_scb[(ind.combn+1):n], epanech_kernel((ind.combn-n):(n-ind.combn)/(n*b))*((ind.combn-n):(n-ind.combn)/(n)),
                                       type="open")[(n-ind.combn+1):(2*(n-ind.combn))]
    dat_fit_ll.right.combn <-   (s2.right.combn*f_tilde.right.combn - s1.right.combn * f_tilde_1.right.combn)/(s2.right.combn * s0.right.combn- s1.right.combn^2)

    dat_fit_ll.left.combn.og<- 2*dat_fit_ll.left.combn - dat_fit_ll.left.combn
    dat_fit_ll.right.combn.og <- 2*dat_fit_ll.right.combn - dat_fit_ll.right.combn

    fit.test.2[,r]<- c(dat_fit_ll.left.combn.og, dat_fit_ll.right.combn.og)
  }

  G <- floor(n*b)
  Sigma<- matrix(0, nrow=d, ncol=d)
  for(i in 1:(n-G+1)){
    Sigma<- Sigma+ colSums(Y[i:(i+G-1),]) %*% t(colSums(Y[i:(i+G-1),]) )
  }
  Sigma<- Sigma/(G*(n-G+1))
  ts.boot2<-c()
  ts.boot1<-c()
  set.seed(1609)
  #st2<-Sys.time()
  for(iter in 1:nboot){
    theta_sup.boot<- c(); theta.cum.boot<-rep(0,right_lim+1-left_lim)
    theta_sup.boot.2<- c(); theta.cum.boot.2<-rep(0,right_lim+1-left_lim)
    Z<- MASS::mvrnorm(n, rep(0,d), Sigma)
    for(r in 1:d){
      theta_val.boot<-c(); theta_val.boot.2<-c()
      dat_scb_1<- Z[,r]+combn_fit.og[,r]

      f_tilde_fft_L <- convolve(dat_scb_1[1:n], epanech_kernel((2*n*b-1):(-n*b)/(n*b)), type="open")[(left_lim):(right_lim)]
      f_tilde1_fft_L <- convolve(dat_scb_1[1:n], epanech_kernel((2*n*b-1):(-n*b)/(n*b))*((2*n*b-1):(-n*b)/(n)), type="open")[(left_lim):(right_lim)]

      dat_fit_ll_L <-   (S2_mat[r, (left_lim:right_lim)-left_lim+1,1]*f_tilde_fft_L -
                           S1_mat[r, (left_lim:right_lim)-left_lim+1,1]*f_tilde1_fft_L)/
        (S2_mat[r, (left_lim:right_lim)-left_lim+1,1]*S0_mat[r, (left_lim:right_lim)-left_lim+1,1]-
           S1_mat[r, (left_lim:right_lim)-left_lim+1,1]^2)

      f_tilde_fft_R <- convolve(dat_scb_1[1:n], epanech_kernel((n*b):(1-2*n*b)/(n*b)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]
      f_tilde1_fft_R <- convolve(dat_scb_1[1:n], epanech_kernel((n*b):(1-2*n*b)/(n*b))*((n*b):(1-2*n*b)/(n)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]

      dat_fit_ll_R <-   (S2_mat[r, (left_lim:right_lim)-left_lim+1,2]*f_tilde_fft_R -
                           S1_mat[r, (left_lim:right_lim)-left_lim+1,2]*f_tilde1_fft_R)/
        (S2_mat[r, (left_lim:right_lim)-left_lim+1,2]*S0_mat[r, (left_lim:right_lim)-left_lim+1,2]-
           S1_mat[r, (left_lim:right_lim)-left_lim+1,2]^2)

      dat_fit_ll_L_og<- 2*dat_fit_ll_L- dat_fit_ll_L
      dat_fit_ll_R_og<- 2*dat_fit_ll_R - dat_fit_ll_R

      theta_val.boot<- abs(dat_fit_ll_L_og-dat_fit_ll_R_og)^2

      dat_scb_2<- Z[,r]+fit.test.2[,r]

      f_tilde_fft.2_L <- convolve(dat_scb_2[1:n], epanech_kernel((2*n*b-1):(-n*b)/(n*b)), type="open")[(left_lim):(right_lim)]
      f_tilde1_fft.2_L <- convolve(dat_scb_2[1:n], epanech_kernel((2*n*b-1):(-n*b)/(n*b))*((2*n*b-1):(-n*b)/(n)), type="open")[(left_lim):(right_lim)]

      dat_fit.2_ll_L <-   (S2_mat[r, (left_lim:right_lim)-left_lim+1,1]*f_tilde_fft.2_L -
                             S1_mat[r, (left_lim:right_lim)-left_lim+1,1]*f_tilde1_fft.2_L)/
        (S2_mat[r, (left_lim:right_lim)-left_lim+1,1]*S0_mat[r, (left_lim:right_lim)-left_lim+1,1]-
           S1_mat[r, (left_lim:right_lim)-left_lim+1,1]^2)

      f_tilde_fft.2_R <- convolve(dat_scb_2[1:n], epanech_kernel((n*b):(1-2*n*b)/(n*b)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]
      f_tilde1_fft.2_R <- convolve(dat_scb_2[1:n], epanech_kernel((n*b):(1-2*n*b)/(n*b))*((n*b):(1-2*n*b)/(n)), type="open")[(left_lim+3*n*b-1):(right_lim+3*n*b-1)]

      dat_fit.2_ll_R <-   (S2_mat[r, (left_lim:right_lim)-left_lim+1,2]*f_tilde_fft.2_R -
                             S1_mat[r, (left_lim:right_lim)-left_lim+1,2]*f_tilde1_fft.2_R)/
        (S2_mat[r, (left_lim:right_lim)-left_lim+1,2]*S0_mat[r, (left_lim:right_lim)-left_lim+1,2]-
           S1_mat[r, (left_lim:right_lim)-left_lim+1,2]^2)

      dat_fit.2_ll_L_og<- 2*dat_fit.2_ll_L - dat_fit.2_ll_L
      dat_fit.2_ll_R_og<- 2*dat_fit.2_ll_R - dat_fit.2_ll_R


      theta_val.boot.2<- abs(dat_fit.2_ll_L_og-dat_fit.2_ll_R_og)^2

      # theta_val.boot<-c()
      # for(i in left_lim:right_lim){
      #   theta_val.boot[i-left_lim +1] <- abs(theta(Z[,r], max(floor(i -  3*n*b+1),1), i,(i- n*b),
      #                                         S0_mat[r, i-left_lim+1,1], S1_mat[r, i-left_lim+1,1], S2_mat[r, i-left_lim+1,1])-
      #                                     theta(Z[,r], i, min(floor(i+3*n*b),n),
      #                                           (i+ n*b),
      #                                           S0_mat[r, i-left_lim+1,2], S1_mat[r, i-left_lim+1,2], S2_mat[r, i-left_lim+1,2] ))
      #
      # }
      theta.cum.boot<- theta_val.boot+theta.cum.boot
      theta_sup.boot[r] <- max(theta_val.boot)

      theta.cum.boot.2<- theta_val.boot.2+theta.cum.boot.2
      theta_sup.boot.2[r] <- max(theta_val.boot.2)
    }
    theta.combn.boot.2<- max(theta.cum.boot.2)
    ts.boot1[iter]<- sum(theta_sup.boot)
    ts.boot2[iter]<- abs(theta.combn.boot.2 - sum(theta_sup.boot.2))
  }
  #et2<-Sys.time()
  dec_1<- (sum(theta_sup)>quantile(ts.boot1, 0.95))
  dec_2<-(ts>quantile(ts.boot2, 0.95))
  return(list(theta_sup=theta_sup, tau_list=tau, tau.combn=tau.combn, ts=ts, theta.combn=theta.combn,
              index.cp=index.cp, dec_1=dec_1, dec_2=dec_2, Sigma=Sigma))
}

