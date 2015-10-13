library(viridis)
CovToCor = function(X){
	# Normalize a covariance matrix into a correlation matrix
    corX=diag(1/sqrt(diag(X))) %*% X %*% diag(1/sqrt(diag(X)));
    return(corX)
}

trace_plot = function(data,main = NULL){
	plot(NA,NA,xlim = c(0,ncol(data)),ylim = range(c(0,data)),main= main)
	for(i in 1:nrow(data)) lines(data[i,],col=i)
}

draw_simulation_diagnostics_new = function(
                                                GxE_state = GxE_state,
                                                sp_num = sp_num,
                                                Posterior    = Posterior,
                                                Lambda       = Lambda,
                                                F            = F,
                                                B_F          = B_F,
                                                F_a          = F_a,
                                                mu           = mu,
                                                B_resid      = B_resid,
                                                E_a2         = E_a2,
                                                E_a          = E_a,
                                                delta        = delta,
                                                B_shape      = B_shape,
                                                F_h2         = F_h2,
                                                prec_B_F     = prec_B_F,
                                                prec_F_resid = prec_F_resid,
                                                resid_Y_prec = resid_Y_prec,
                                                resid_h2     = resid_h2,
                                                E_a_prec     = E_a_prec,
                                                E_a2_prec    = E_a2_prec
                ){

    devices = dev.list()
    while(length(devices) < 5){
    	quartz()
    	devices = dev.list()
    }
    dev.set(devices[1])
    par(mfrow=c(3,3))

    p = nrow(Lambda)


    sim_data = GxE_state$run_parameters$setup
    actual_E_Lambda = sim_data$Lambda

    cors = abs(cor(Lambda,actual_E_Lambda))
    cors=rbind(cors,apply(cors,2,max))
    cors = cbind(cors,apply(cors,1,max))
    image(cors[,ncol(cors):1],zlim=c(0,1))

    B = B_resid + B_F %*% t(Lambda)
    B_act = sim_data$B + sim_data$B_F %*% t(sim_data$Lambda)
    plot(B,B_act);abline(0,1)
    plot(B_F %*% t(Lambda),B_act);abline(0,1)
    plot(c(B_resid,B_F),pch=21,bg=c(rep(1,length(B_resid)),rep(2,length(B_F))))
    plot(c(B_F %*% t(Lambda)),pch=21,bg=c(rep(1,length(B_resid)),rep(2,length(B_F))))

    plot(mu,sim_data$mu)

    plot(F_h2,type='l',ylim=c(0,1))
    lines(1-F_h2,col=2)

    var_B_F = 1/prec_B_F
    var_F_a = 1/prec_F_resid * F_h2
    var_F_r = 1/prec_F_resid * (1-F_h2)
    image(t(cbind(var_B_F,var_F_a,var_F_r)))

    if(sp_num>1){

    	dev.set(devices[2])
        # Figure of trace plots of the largest elements of each column of the factor loading matrix
        f2_row=4;
        f2_col=4;
        par(mfrow=c(f2_row,f2_col))
        for(k in 0:(min(2*f2_row,nrow(Posterior$Lambda)/p)-1)) {
            o = order(-abs(rowMeans(Posterior$Lambda[(1:p)+k*p,max(1,sp_num-100):sp_num])))
            traces = Posterior$Lambda[o[1:5]+k*p,1:sp_num]
            trace_plot(traces)
            
        }

    	dev.set(devices[3])
    	# Figure of trace plots and histograms of the factor heritablities
        f4_row=4;
        f4_col=4;
        par(mfrow=c(f4_row,f4_col))
        for(k in 1:min(2*f4_row,nrow(Posterior$F_h2))) {
            h2s = Posterior$F_h2[k,1:sp_num]
            if(sum(h2s[!is.na(h2s)])==0) {
                next
            }
            plot(h2s,type='l')
            hist(h2s,breaks=100,xlim=c(-0.1,1))
        }

        dev.set(devices[4])
        par(mfrow=c(1,2))
        B_F_L_mean = array(0,dim = dim(B_act))
        B_F_mean = 0*matrix(Posterior$B_F[,1],nr = dim(B)[1])
        B_resid_mean = 0*B_resid
        for(i in 1:sp_num){
            Lambda_i = matrix(Posterior$Lambda[,i],nr = p)
            B_F_i = matrix(Posterior$B_F[,i],nr = dim(B)[1])
            B_F_L_mean = B_F_L_mean + (B_F_i %*% t(Lambda_i)) / sp_num
            B_resid_i = matrix(Posterior$B_resid[,i],nr = dim(B_resid)[1])
            B_F_mean = B_F_mean + B_F_i/sp_num
            B_resid_mean = B_resid_mean + B_resid_i/sp_num
        }
        plot(B_F_L_mean,B_act);abline(0,1)
        plot(c(B_resid_mean,B_F_mean),pch=21,bg=c(rep(1,length(B_resid_mean)),rep(2,length(B_F_mean))))

        dev.set(devices[5])
        k = nrow(Posterior$Lambda)/p
        Lambda = matrix(0,p,k)
        for(i in 1:k){
            j = p*(i-1)+1:p
            Lambda = Lambda + matrix(Posterior$Lambda,p,k)/I
        }
        Lambda_thresh = apply(Lambda,2,function(x) {x[abs(x)<quantile(abs(x),.75)]=0;x})
        Lambda = Lambda[order(Lambda_thresh[,1],Lambda_thresh[,2],Lambda_thresh[,3],Lambda_thresh[,4],Lambda_thresh[,5],Lambda_thresh[,6]),]
        image(t((Lambda)),col=viridis(250))


     #    k = nrow(Posterior$Lambda)/p;
     #    h2s = Posterior$F_h2[,1:sp_num]
     #    G_Lambdas = array(0,dim = dim(Posterior$Lambda))
     #    Lambda_est = matrix(0,p,k)
     #    G_est = E_est = matrix(0,p,p)
     #    for(j in 1:sp_num) {
     #        Lj = matrix(Posterior$Lambda[,j],p,k);
     #        h2j = Posterior$F_h2[,j];
     #        G_Lj = Lj %*%  diag(sqrt(h2j));
     #        G_Lambdas[,j] = c(G_Lj)
     #        Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
     #        G_est = G_est + Gj/sp_num
            
     #        E_Lj = Lj  %*% diag(1-h2j) %*%  t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
     #        E_est = E_est + E_Lj/sp_num;
     #        Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
     #    }
     #    G_Lambda = matrix(rowMeans(G_Lambdas),p,k)

     #    dev.set(devices[4])
     #    par(mfrow=c(3,3))
     #    cors = abs(cor(G_Lambda,actual_E_Lambda))
    	# cors=rbind(cors,apply(cors,2,max))
    	# cors = cbind(cors,apply(cors,1,max))
    	# image(cors[,ncol(cors):1],zlim=c(0,1))

    	# image(CovToCor(G_est),zlim=c(-1,1))
    	# image(CovToCor(E_est),zlim=c(-1,1))

    	# clim=max(0.1,min(.75,max(max(abs(G_Lambda)))))
    	# clims=c(-clim,clim)
    	# image(t(G_Lambda),zlim=clims)

    	# image(t(actual_G_Lambda),zlim=clims)

    	# plot(c(G_est),c(G_act));abline(0,1)
     #    plot(c(E_est),c(E_act));abline(0,1)

     #    plot(diag(G_est)/(diag(G_est)+diag(E_est)),h2s_act,xlim=c(0,1),ylim=c(0,1))
     #    abline(0,1)

    	# F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
    	# E_h2 = 1-F_h2
    	# plot(F_h2,type='l',ylim=c(0,1))
    	# lines(E_h2,col=2)
    }
}

draw_results_diagnostics_new = function(
                                                GxE_state = GxE_state,
                                                sp_num = sp_num,
                                                Posterior    = Posterior,
                                                Lambda       = Lambda,
                                                F            = F,
                                                B_F          = B_F,
                                                F_a          = F_a,
                                                mu           = mu,
                                                B_resid      = B_resid,
                                                E_a2         = E_a2,
                                                E_a          = E_a,
                                                delta        = delta,
                                                B_shape      = B_shape,
                                                F_h2         = F_h2,
                                                prec_B_F     = prec_B_F,
                                                prec_F_resid = prec_F_resid,
                                                resid_Y_prec = resid_Y_prec,
                                                resid_h2     = resid_h2,
                                                E_a_prec     = E_a_prec,
                                                E_a2_prec    = E_a2_prec
                ){

    devices = dev.list()
    while(length(devices) < 4){
        quartz()
        devices = dev.list()
    }
    dev.set(devices[1])
    par(mfrow=c(3,3))

    p = nrow(Lambda)

    plot(F_h2,type='l',ylim=c(0,1))
    lines(1-F_h2,col=2)

    var_B_F = 1/prec_B_F
    var_F_a = 1/prec_F_resid * F_h2
    var_F_r = 1/prec_F_resid * (1-F_h2)
    image(t(cbind(var_B_F,var_F_a,var_F_r)))

    if(sp_num>1){

        dev.set(devices[2])
        # Figure of trace plots of the largest elements of each column of the factor loading matrix
        f2_row=4;
        f2_col=4;
        par(mfrow=c(f2_row,f2_col))
        for(k in 0:(min(2*f2_row,nrow(Posterior$Lambda)/p)-1)) {
            o = order(-abs(rowMeans(Posterior$Lambda[(1:p)+k*p,max(1,sp_num-100):sp_num])))
            traces = Posterior$Lambda[o[1:5]+k*p,1:sp_num]
            trace_plot(traces)
            
        }

        dev.set(devices[3])
        # Figure of trace plots and histograms of the factor heritablities
        f4_row=4;
        f4_col=4;
        par(mfrow=c(f4_row,f4_col))
        for(k in 1:min(2*f4_row,nrow(Posterior$F_h2))) {
            h2s = Posterior$F_h2[k,1:sp_num]
            if(sum(h2s[!is.na(h2s)])==0) {
                next
            }
            plot(h2s,type='l',ylim=c(-0.1,1))
            hist(h2s,breaks=100,xlim=c(-0.1,1))
        }

        dev.set(devices[4])
        # Figure of trace plots and histograms of the factor b
        f4_row=4;
        f4_col=4;
        par(mfrow=c(f4_row,f4_col))
        for(k in 1:min(2*f4_row,nrow(Posterior$B_F))) {
            bs = Posterior$B_F[k,1:sp_num]
            if(sum(bs[!is.na(bs)])==0) {
                next
            }
            plot(bs,type='l',ylim=range(c(0,bs)));abline(h=0)
            hist(bs,breaks=100,xlim=range(c(0,bs)));abline(v=0)
        }


        dev.set(devices[5])
        k = nrow(Posterior$Lambda)/p
        Lambda = matrix(0,p,k)
        for(i in 1:k){
            j = p*(i-1)+1:p
            Lambda = Lambda + matrix(Posterior$Lambda,p,k)/sp_num
        }
        Lambda_thresh = apply(Lambda,2,function(x) {x[abs(x)<quantile(abs(x),.75)]=0;x})
        p_order = order(Lambda_thresh[,1],Lambda_thresh[,2],Lambda_thresh[,3],Lambda_thresh[,4],Lambda_thresh[,5],Lambda_thresh[,6])
        Lambda = Lambda[p_order,]
        image(t((Lambda)),col=viridis(250))



        # par(mfrow=c(1,2))
        # B_F_L_mean = array(0,dim = dim(B_act))
        # B_F_mean = 0*matrix(Posterior$B_F[,1],nr = dim(B)[1])
        # B_resid_mean = 0*B_resid
        # for(i in 1:sp_num){
        #     Lambda_i = matrix(Posterior$Lambda[,i],nr = p)
        #     B_F_i = matrix(Posterior$B_F[,i],nr = dim(B)[1])
        #     B_F_L_mean = B_F_L_mean + (B_F_i %*% t(Lambda_i)) / sp_num
        #     B_resid_i = matrix(Posterior$B_resid[,i],nr = dim(B_resid)[1])
        #     B_F_mean = B_F_mean + B_F_i/sp_num
        #     B_resid_mean = B_resid_mean + B_resid_i/sp_num
        # }
        # plot(B_F_L_mean,B_act);abline(0,1)
        # plot(c(B_resid_mean,B_F_mean),pch=21,bg=c(rep(1,length(B_resid_mean)),rep(2,length(B_F_mean))))

        k = nrow(Posterior$Lambda)/p;
        h2s = Posterior$F_h2[,1:sp_num]
        G_Lambdas = array(0,dim = dim(Posterior$Lambda))
        Lambda_est = matrix(0,p,k)
        G_est = E_est = P_est = matrix(0,p,p)
        for(j in 1:sp_num) {
            Lj = matrix(Posterior$Lambda[,j],p,k);
            s2aj = Posterior$F_h2[,j]/Posterior$prec_F_resid[,j];
            s2ej = (1-Posterior$F_h2[,j])/Posterior$prec_F_resid[,j];
            s2pj = 1/Posterior$prec_F_resid[,j];
            s2aj[is.na(s2aj)] = 0
            s2ej[is.na(s2ej)] = 0
            s2pj[is.na(s2pj)] = 0
            s2aj[s2aj==Inf] = 0
            s2ej[s2ej==Inf] = 0
            s2pj[s2pj==Inf] = 0
            G_Lj = Lj %*%  diag(sqrt(s2aj));
            G_Lambdas[,j] = c(G_Lj)
            Gj = G_Lj %*%  t(G_Lj) + diag(Posterior$resid_h2[,j]/Posterior$resid_Y_prec[,j])
            G_est = G_est + Gj/sp_num
            
            E_Lj = Lj  %*% diag(s2ej) %*%  t(Lj) + diag((1-Posterior$resid_h2[,j])/Posterior$resid_Y_prec[,j])
            E_est = E_est + E_Lj/sp_num;
            Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
            # Pj = Lj  %*% diag(s2pj)%*% t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
            # P_est = P_est + Pj/sp_num
        }
        G_Lambda = matrix(rowMeans(G_Lambdas),p,k)
        image(CovToCor(G_est[p_order,p_order]),zlim=c(-1,1),col=viridis(250))
        image(CovToCor(E_est[p_order,p_order]),zlim=c(-1,1),col=viridis(250))
        # image(CovToCor(P_est[p_order,p_order]),zlim=c(-1,1),col=viridis(250))
        # plot(P_est,cov(GxE_state$data$Y_full));abline(0,1)
        plot(G_est + E_est,cov(GxE_state$data$Y_full));abline(0,1)
        # plot(G_est + E_est,P_est);abline(0,1)

     #    dev.set(devices[4])
     #    par(mfrow=c(3,3))
     #    cors = abs(cor(G_Lambda,actual_E_Lambda))
        # cors=rbind(cors,apply(cors,2,max))
        # cors = cbind(cors,apply(cors,1,max))
        # image(cors[,ncol(cors):1],zlim=c(0,1))

        # image(CovToCor(G_est),zlim=c(-1,1))
        # image(CovToCor(E_est),zlim=c(-1,1))

        # clim=max(0.1,min(.75,max(max(abs(G_Lambda)))))
        # clims=c(-clim,clim)
        # image(t(G_Lambda),zlim=clims)

        # image(t(actual_G_Lambda),zlim=clims)

        # plot(c(G_est),c(G_act));abline(0,1)
     #    plot(c(E_est),c(E_act));abline(0,1)

     #    plot(diag(G_est)/(diag(G_est)+diag(E_est)),h2s_act,xlim=c(0,1),ylim=c(0,1))
     #    abline(0,1)

        # F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
        # E_h2 = 1-F_h2
        # plot(F_h2,type='l',ylim=c(0,1))
        # lines(E_h2,col=2)
    }
}
