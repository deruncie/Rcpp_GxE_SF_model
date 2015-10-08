

diagon2 = function(x,y) {
    sy = svd(y)
    K = sy$u %*% diag(sqrt(sy$d))
    Kinv = solve(K)
    G = Kinv %*% x %*% t(Kinv)
    eG = svd(G)
    S = diag(1/sqrt(eG$d)) %*% t(eG$u) %*% Kinv
    return(S)
}


GSVD_2 = function(A,B){
	library(pracma)
	# There is almost no advantage to going to RcppArmadillo for this function
	# UCt(X) = A
	# VSt(X) = B
	# [U,d,V] = svd(A*inv(B)): UCt(X) * inv(t(X)) * inv(S) * t(V) = A*inv(B)
	# solve: Ct(X) = t(U)*A, St(X) = t(V)*B. Set S = diag(1,n), C = C*inv(S). Then Normalize C^2 + S^2 + Eye(n), and correspondingly normalize columns of X.
	
	dimA = dim(A)
	dimB = dim(B)
	AiB = A %*% pinv(B)
	svd1 = svd(AiB)

	recover()

	U = svd1$u
	V = svd1$v

	d = svd1$d

	x1 = t(B) %*% V
	c1 = d
	s1 = rep(1,length(c1))


	# X = sweep(x1,2,norm_factor,'*')

	if(dimB[1]<dimB[2]){
		r = diff(dimB)
		Aresid = A - U %*% diag(c1) %*% t(x1)
		svd2 = svd(Aresid)
		U2 = cbind(svd2$u[,1:r],U)
		XA =svd2$v[,1:r]
		A - U2 %*% diag(c(svd2$d[1:r],c1)) %*% t(cbind(XA,x1))
		c1 = c(rep(0,r),c1)
		X = cbind(XA,x1)
	}
	if(dimA[1]<dimA[2]){
		r = diff(dimA)
		Bresid = B - V %*% diag(s1) %*% t(x1)
		svd2 = svd(Bresid)
		V2 = cbind(V,svd2$u[,1:r])
		XB = svd2$v[,1:r]
		B - V2 %*% diag(c(s1,svd2$d[1:r])) %*% t(cbind(x1,XB))
		s1 = c(s1,rep(0,r))
		X = cbind(X,XB)
	}
		


	norm_factor = sqrt(c1^2 + s1^2)
	c2 = c1 / norm_factor
	s2 = s1 / norm_factor

	return(list(U=U,V=V,X=X,C=C,S=S))
}

A = matrix(rnorm(4*6),4);A = t(A) %*% A
a = cholcov(A)
B = matrix(rnorm(3*6),3);B = t(B) %*% B
b = cholcov(B)
GSVD_2(a,b)
GSVD(a,b)

library(R.matlab)
Matlab$startServer("/Applications/MATLAB_R2014a.app/bin/matlab")
matlab <- Matlab()
isOpen <- open(matlab)
setVariable(matlab, a=a)
setVariable(matlab, b=b)
evaluate(matlab,'[U,V,X,C,S]=gsvd(a,b)')
result = list()
result$U <- getVariable(matlab, "U")$U
result$V <- getVariable(matlab, "V")$V
result$X <- getVariable(matlab, "X")$X
result$C <- getVariable(matlab, "C")$C
result$S <- getVariable(matlab, "S")$S
test_GSVD(a,b,result)
test_GSVD = function(A,B,result){
	print(result$U %*% result$C %*% t(result$X) - A)
	print(result$V %*% result$S %*% t(result$X) - B)
}

library(MASS)
test_GSVD2 = function(A,B,result){
	s1 = diag(result$C)^2
	s2 = diag(result$S)^2
	if(dim(result$X)[2] < dim(result$X)[1]){
		null = Null(result$X)
		result$X = cbind(result$X,null)
		s1 = c(s1,Inf)
		s2 = c(s2,Inf)
	}
	U = t(solve(result$X))

	return(U %*% diag(1/(3*s1+5*s2)) %*% t(U) - solve(3*t(A) %*% A + 5*t(B) %*% B))
}

test_GSVD3 = function(A,B,result){
	U = result$X
	s1 = diag(result$C)^2
	s2 = diag(result$S)^2

	return(U %*% diag((3*s1+5*s2)) %*% t(U) - (3*t(A) %*% A + 5*t(B) %*% B))
}



GSVD_PEIP = function (A, B) 
{
    if (!is.matrix(A)) 
        stop("Argument A should be a matrix")
    if (!is.matrix(B)) 
        stop("Argument B should be a matrix")
    dimA = dim(A)
    dimB = dim(B)
    if (dimA[1] == 0) 
        stop("Matrix A has zero rows/columns")
    if (dimB[1] == 0) 
        stop("Matrix B has zero rows/columns")
    if (!all(is.finite(A))) 
        stop("Matrix A may not contain infinite/NaN/NA")
    if (!all(is.finite(B))) 
        stop("Matrix B may not contain infinite/NaN/NA")
    ncolA = dimA[2]
    nrowA = dimA[1]
    ncolB = dimB[2]
    nrowB = dimB[1]
    lwork = max(c(3 * ncolA, nrowA, nrowB)) + ncolA
    ju = 2
    jv = 2
    jq = 2
    z <- .Fortran("zdggsvd", as.integer(ju), as.integer(jv), 
        as.integer(jq), as.integer(nrowA), as.integer(ncolA), 
        as.integer(nrowB), integer(1), integer(1), as.double(A), 
        as.integer(nrowA), as.double(B), as.integer(nrowB), double(ncolA), 
        double(ncolA), double(nrowA * nrowA), as.integer(nrowA), 
        double(nrowB * nrowB), as.integer(nrowB), double(ncolA * 
            ncolA), as.integer(ncolA), double(lwork), integer(ncolA), 
        integer(1), dup = FALSE, PACKAGE = "PEIP")
    K = z[7][[1]]
    L = z[8][[1]]
    U = z[15][[1]]
    V = z[17][[1]]
    Q = z[19][[1]]
    ALPHA = z[13][[1]]
    BETA = z[14][[1]]
    R = matrix(z[9][[1]], ncol(A), nrow = nrow(A), byrow = FALSE)
    U = matrix(U, ncol = nrow(A), nrow = nrow(A), byrow = FALSE)
    V = matrix(V, ncol = nrow(B), nrow = nrow(B), byrow = FALSE)
    Q = matrix(Q, ncol = ncol(A), nrow = ncol(A), byrow = FALSE)
    D1 = mat.or.vec(nrow(A), K + L)
    D2 = mat.or.vec(nrow(B), K + L)
    oR = mat.or.vec((K + L), ncol(A))
    # recover()
    if (K > 0) {
        if (K == 1) {
            D1[1:K, 1:K] = rep(1, K)
        }
        else {
            diag(D1[1:K, 1:K]) = rep(1, K)
        }
        diag(D1[,1:nrow(D1)]) = ALPHA[1:nrow(D1)]
        # diag(D1[(K + 1):(K + L), (K + 1):(K + L)]) = ALPHA[(K + 1):(K + L)]
        diag(D2[1:L, (K + 1):(K + L)]) = BETA[(K + 1):(K + L)]
    }
    if (K == 0) {
        diag(D1[(K + 1):(K + L), (K + 1):(K + L)]) = ALPHA[(K + 
            1):(K + L)]
        diag(D2[1:L, (K + 1):(K + L)]) = BETA[(K + 1):(K + L)]
    }
    Ci = ALPHA[(K + 1):(K + L)]
    S = BETA[(K + 1):(K + L)]
    oR = R
    # oR[(1):(K + L), (ncol(A) - K - L + 1):(ncol(A))] = R[, (ncol(A) - K - L + 1):(ncol(A))]
    X = t(oR %*% t(Q))
    return(list(U = U, V = V, X = X, C = D1, S = D2))
}