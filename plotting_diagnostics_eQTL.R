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
                                                eQTL_state   = eQTL_state,
                                                sp_num       = sp_num,
                                                Posterior    = Posterior,
                                                Lambda       = Lambda,
                                                F            = F,
                                                B_F          = B_F,
                                                E_a2_F       = E_a2_F,
                                                F_a          = F_a,
                                                mu           = mu,
                                                B_resid      = B_resid,
                                                cis_effects  = cis_effects,
                                                E_a2         = E_a2,
                                                E_a          = E_a,
                                                delta        = delta,
                                                B_shape      = B_shape,
                                                F_h2         = F_h2,
                                                prec_B_F     = prec_B_F,
                                                prec_E_a2_F  = prec_E_a2_F,
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


    sim_data = eQTL_state$run_parameters$setup
    actual_E_Lambda = sim_data$Lambda

    cors = abs(cor(Lambda,actual_E_Lambda))
    cors=rbind(cors,apply(cors,2,max))
    cors = cbind(cors,apply(cors,1,max))
    image(cors[,ncol(cors):1],zlim=c(0,1),col=viridis(250))

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
    image(t(cbind(var_B_F,var_F_a,var_F_r)),col = viridis(250))

    plot(cis_effects,sim_data$cis_effects,xlim=c(-4,4),ylim=c(-4,4));abline(0,1)

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
        par(mfrow=c(3,3))
        k = nrow(Posterior$Lambda)/p
        Lambda = matrix(0,p,k)
        for(i in 1:k){
            j = p*(i-1)+1:p
            Lambda = Lambda + matrix(Posterior$Lambda,p,k)/sp_num
        }
        Lambda_thresh = apply(Lambda,2,function(x) {x[abs(x)<quantile(abs(x),.75)]=0;x})
        p_order = order(Lambda_thresh[,1],Lambda_thresh[,2],Lambda_thresh[,3],Lambda_thresh[,4],Lambda_thresh[,5],Lambda_thresh[,6])
        Lambda = Lambda[p_order,]
        image(t((Lambda)),col=viridis(250),zlim = c(-1,1)*max(abs(Lambda)))


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
        plot(G_est + E_est,cov(eQTL_state$data$Y_full));abline(0,1)

        dev.set(devices[5])
        boxplot(t(Posterior$cis_effects),ylim=c(-10,10))
        plot(rowMeans(Posterior$cis_effects),sim_data$cis_effects,xlim=c(-4,4),ylim=c(-4,4));abline(0,1)

        F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
        E_h2 = 1-F_h2
        plot(F_h2,type='l',ylim=c(0,1))
        lines(E_h2,col=2)
    }
}

draw_results_diagnostics_new = function(
                                                eQTL_state   = eQTL_state,
                                                sp_num       = sp_num,
                                                Posterior    = Posterior,
                                                Lambda       = Lambda,
                                                F            = F,
                                                B_F          = B_F,
                                                E_a2_F       = E_a2_F,
                                                F_a          = F_a,
                                                mu           = mu,
                                                B_resid      = B_resid,
                                                cis_effects  = cis_effects,
                                                E_a2         = E_a2,
                                                E_a          = E_a,
                                                delta        = delta,
                                                B_shape      = B_shape,
                                                F_h2         = F_h2,
                                                prec_B_F     = prec_B_F,
                                                prec_E_a2_F  = prec_E_a2_F,
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
    par(mfrow=c(2,2))

    p = nrow(Posterior$mu)
    sp_num = ncol(Posterior$mu)

    plot(F_h2,type='l',ylim=c(0,1))
    lines(1-F_h2,col=2)

    var_B_F = 1/prec_B_F
    var_E_a2_F = 1/prec_E_a2_F
    var_F_a = 1/prec_F_resid * F_h2
    var_F_r = 1/prec_F_resid * (1-F_h2)
    image(t(cbind(var_B_F,var_E_a2_F,var_F_a,var_F_r)),col=viridis(250))

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

        if(nrow(Posterior$E_a2_F)>0){
            while(length(devices) < 6){
                quartz()
                devices = dev.list()
            }
            dev.set(devices[6])
            # Figure of trace plots and histograms of the factor heritablities
            f4_row=4;
            f4_col=4;
            par(mfrow=c(f4_row,f4_col))
            for(k in 1:min(f4_col*f4_row,nrow(Posterior$prec_E_a2_F))) {
                s2a2 = 1/Posterior$prec_E_a2_F[k,1:sp_num]
                if(sum(h2s[!is.na(h2s)])==0) {
                    next
                }
                plot(s2a2,type='l',ylim=c(-0.1,1));abline(h=0)
            }
            # par(mfrow=c(1,1))
            # k = nrow(Posterior$Lambda)/p
            # j = 6
            # r = eQTL_state$run_variables$r
            # F_a_s = t(Posterior$F_a[(j-1)*r + 1:r,])
            # F_a_s = F_a_s[,order(colMeans(F_a_s))]
            # boxplot(F_a_s,main = j);abline(h=0)
            # n = eQTL_state$run_variables$n
            # F_s = t(Posterior$F[(j-1)*n + 1:n,])
            # F_s = F_s[,order(colMeans(F_s))]
            # boxplot(F_s,main = j);abline(h=0)
        }


        dev.set(devices[5])
        par(mfrow=c(2,2))
        k = nrow(Posterior$Lambda)/p
        Lambda = matrix(0,p,k)
        for(i in 1:k){
            j = p*(i-1)+1:p
            Lambda = Lambda + matrix(Posterior$Lambda,p,k)/sp_num
        }
        Lambda_thresh = apply(Lambda,2,function(x) {x[abs(x)<quantile(abs(x),.75)]=0;x})
        p_order = order(Lambda_thresh[,1],Lambda_thresh[,2],Lambda_thresh[,3],Lambda_thresh[,4],Lambda_thresh[,5],Lambda_thresh[,6])
        Lambda = Lambda[p_order,]
        image(t((Lambda)),col=viridis(250),zlim = c(-1,1)*max(abs(Lambda)))


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
        plot(G_est + E_est,cov(eQTL_state$data$Y_full));abline(0,1)

        dev.set(devices[1])
        par(mfrow=c(2,2))
        # recover()
        boxplot(t(1/Posterior$prec_E_a2_F))
        boxplot(t(1/Posterior$prec_F_resid * Posterior$F_h2))
        boxplot(t(1/Posterior$prec_F_resid * (1-Posterior$F_h2)))


        genes_with_cisG = which(colSums(eQTL_state$data_matrices$cisGenotypes!=0)>0)
        boxplot(t(Posterior$cis_effects[genes_with_cisG,]));abline(h=0)

        boxplot(t(Posterior$B_F))

        boxplot(1/t(apply(Posterior$delta,2,function(x) cumprod(x))))
    }
}



posterior_plots = function(Posterior){

    devices = dev.list()
    while(length(devices) < 5){
        quartz()
        devices = dev.list()
    }
    dev.set(devices[1])
    par(mfrow=c(2,2))

    p = nrow(Posterior$mu)
    sp_num = ncol(Posterior$mu)

    if(sp_num>1){

        lambda_vars = rowMeans(1/apply(Posterior$delta,2,cumprod))
        lambda_vars[lambda_vars==Inf]=NA
        plot(lambda_vars)
        factor_order = order(-lambda_vars)
        major_factors = na.omit(factor_order[lambda_vars[factor_order] > 0.01])

        dev.set(devices[2])
        # Figure of trace plots of the largest elements of each column of the factor loading matrix
        f2_row=4;
        f2_col=4;
        par(mfrow=c(f2_row,f2_col))
        for(k in major_factors) {
            o = order(-abs(rowMeans(Posterior$Lambda[(k-1)*p+(1:p),max(1,sp_num-100):sp_num])))
            traces = Posterior$Lambda[(k-1)*p+o[1:5],1:sp_num]
            trace_plot(traces)            
        }

        dev.set(devices[3])
        # Figure of trace plots and histograms of the factor heritablities
        f4_row=4;
        f4_col=4;
        par(mfrow=c(f4_row,f4_col))
        for(k in major_factors) {
            h2s = Posterior$F_h2[k,1:sp_num]
            if(sum(h2s[!is.na(h2s)])==0) {
                next
            }
            plot(h2s,type='l',ylim=c(-0.1,1),main=k)
            hist(h2s,breaks=100,xlim=c(-0.1,1))
        }

        dev.set(devices[4])
        # Figure of trace plots and histograms of the factor b
        f4_row=4;
        f4_col=4;
        par(mfrow=c(f4_row,f4_col))
        for(k in major_factors) {
            bs = Posterior$B_F[k,1:sp_num]
            if(sum(bs[!is.na(bs)])==0) {
                next
            }
            plot(bs,type='l',ylim=range(c(0,bs)),main=k);abline(h=0)
            hist(bs,breaks=100,xlim=range(c(0,bs)));abline(v=0)
        }

        if(nrow(Posterior$E_a2_F)>0){
            while(length(devices) < 6){
                quartz()
                devices = dev.list()
            }
            dev.set(devices[6])
            # Figure of trace plots and histograms of the factor heritablities
            f4_row=4;
            f4_col=4;
            par(mfrow=c(f4_row,f4_col))
            for(k in major_factors) {
                s2a2 = 1/Posterior$prec_E_a2_F[k,1:sp_num]
                if(sum(h2s[!is.na(h2s)])==0) {
                    next
                }
                plot(s2a2,type='l',ylim=c(-0.1,1),main=k);abline(h=0)
            }
        }
            # par(mfrow=c(1,1))
            # k = nrow(Posterior$Lambda)/p
            # j = 6
            # r = eQTL_state$run_variables$r
            # F_a_s = t(Posterior$F_a[(j-1)*r + 1:r,])
            # F_a_s = F_a_s[,order(colMeans(F_a_s))]
            # boxplot(F_a_s,main = j);abline(h=0)
            # n = eQTL_state$run_variables$n
            # F_s = t(Posterior$F[(j-1)*n + 1:n,])
            # F_s = F_s[,order(colMeans(F_s))]
            # boxplot(F_s,main = j);abline(h=0)
        


        dev.set(devices[5])
        par(mfrow=c(2,2))
        k = nrow(Posterior$Lambda)/p
        Lambda = apply(array(Posterior$Lambda,dim=c(p,k,sp_num)),c(1,2),mean)
        Lambda_thresh = apply(Lambda,2,function(x) {x[abs(x)<quantile(abs(x),.75)]=0;x})
        p_order = order(Lambda_thresh[,major_factors[1]],
                        Lambda_thresh[,major_factors[2]],
                        Lambda_thresh[,major_factors[3]],
                        Lambda_thresh[,major_factors[4]],
                        Lambda_thresh[,major_factors[5]],
                        Lambda_thresh[,major_factors[6]],
                        Lambda_thresh[,major_factors[7]],
                        Lambda_thresh[,major_factors[8]]
                        )
        Lambda = Lambda[p_order,major_factors]
        image(t(Lambda),col=viridis(250),zlim = c(-1,1)*max(abs(Lambda)))

        ds = apply(Lambda,2,function(x) density(abs(x)))
        ylim = c(0,max(sapply(ds,function(x)x$y)))
        plot(NA,NA,xlim = c(0,1),ylim=ylim)
        for(i in 1:length(ds)) lines(ds[[i]],col=viridis(10)[i])


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
        plot(G_est + E_est,cov(eQTL_state$data$Y_full));abline(0,1)

        dev.set(devices[1])
        par(mfrow=c(2,2))
        # recover()
        boxplot(t(1/Posterior$prec_E_a2_F))
        boxplot(t(1/Posterior$prec_F_resid * Posterior$F_h2))
        boxplot(t(1/Posterior$prec_F_resid * (1-Posterior$F_h2)))

        par(mfcol=c(3,2))
        library(mcmcplots)
        image(t(Lambda),col=viridis(250),zlim = c(-1,1)*max(abs(Lambda)))
        var_data = t(1/Posterior$prec_F_resid * Posterior$F_h2)[,major_factors];colnames(var_data) = paste('F',1:ncol(var_data),sep='_')
        p = caterplot(var_data,horiz=F,reorder=F,main = 's^2_a')
        # var_data = t(1/Posterior$prec_E_a2_F)[,major_factors];colnames(var_data) = paste('F',1:ncol(var_data),sep='_')
        # p = caterplot(var_data,horiz=F,reorder=F)
        # var_data = t(1/Posterior$prec_F_resid * (1-Posterior$F_h2))[,major_factors];colnames(var_data) = paste('F',1:ncol(var_data),sep='_')
        # p = caterplot(var_data,horiz=F,reorder=F)
        var_data = apply(Posterior$B_F,1,function(x) x*sign(mean(x)))[,major_factors];colnames(var_data) = paste('F',1:ncol(var_data),sep='_')
        p = caterplot(var_data,horiz=F,reorder=F);abline(h=0)

        genes_with_cisG = which(colSums(eQTL_state$data_matrices$cisGenotypes!=0)>0)
        boxplot(t(Posterior$cis_effects[genes_with_cisG,]));abline(h=0)

        boxplot(t(Posterior$B_F))

        boxplot(1/t(apply(Posterior$delta,2,function(x) cumprod(x))))
    }
}


