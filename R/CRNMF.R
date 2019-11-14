

#' @title construct contigency matrix
#' @description This fuction is used to construct contigency matrix.
#' @param mem1 first label
#' @param mem2 second label
#' @return
#' \item{cont}{contigency matrix}
#' @keywords contigency matrix
#' @examples 
#' \dontrun{
#' C<-contigency(mem1,mem2)
#' }
#' @export

contigency<-function(mem1,mem2){
    cont<-matrix(0,nrow=max(mem1),ncol=max(mem2))
    for(i in 1:length(mem1)){
        cont[mem1[i],mem2[i]]<-cont[mem1[i],mem2[i]]+1
    }
    return(cont)
}



#' @title compute Rand index
#' @description This fuction is used to compute Rand index.
#' @param c1 first label
#' @param c2 second label
#' @return
#' \item{AR}{adjusted Rand index}
#' \item{RI}{Rand index}
#' @keywords Rand index
#' @examples 
#' \dontrun{
#' res<-rand.idx(c1,c2)
#' }
#' @export

rand.idx<-function(c1,c2){
    C<-contigency(c1,c2)
    n<-sum(sum(C))
    nis<-sum(rowSums(C)^2)
    njs<-sum(colSums(C)^2)
    t1<-n*(n-1)/2
    t2<-sum(sum(C^2))
    t3<-0.5*(nis+njs)
    nc<-(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1))
    A<-t1+t2-t3
    D<- -t2+t3
    if(t1==nc){
        AR<-0
    }
    else{
        AR<-(A-nc)/(t1-nc)
    }
    RI<-A/t1		
    MI<-D/t1		
    HI<-(A-D)/t1
    return(list(AR=AR,RI=RI))
}



#' @title compute F-score, precision and recall
#' @description This fuction is used to compute F-score, precision and recall.
#' @param G first label
#' @param H second label
#' @return
#' \item{f}{F-score}
#' \item{p}{precision}
#' \item{r}{recall}
#' @keywords F-score, precision, recall
#' @examples 
#' \dontrun{
#' res<-compute.f(c1,c2)
#' }
#' @export

compute.f<-function(G,H){
    N<-length(G)
    numG<-0
    numH<-0
    numI<-0
    for(n in 1:(N-1)){
        Gn <- (G[(n+1):N]==G[n])
        Hn <- (H[(n+1):N]==H[n])
        numG <- numG + sum(Gn)
        numH <- numH + sum(Hn)
        numI <- numI + sum(Gn*Hn)
    }
    p<-1
    r<-1
    f<-1
    if(numH > 0){
        p<-numI / numH
    }
    if(numG > 0){
        r<-numI / numG
    }
    if((p+r) == 0){
        f<-0
    }
    else{
        f<-2 * p * r / (p + r)
    }
    return(list(f=f,p=p,r=r))   
}



#' @title compute normalized mutual information
#' @description This fuction is used to compute normalized mutual information.
#' @param label first label
#' @param result second label
#' @return
#' \item{v}{normalized mutual information}
#' @keywords normalized mutual information
#' @examples 
#' \dontrun{
#' res<-nmi(c1,c2)
#' }
#' @export

nmi<-function(label,result){
    n<-length(label)
    
    label_unique<-sort(unique(label))
    result_unique<-sort(unique(result))
    c<-length(label_unique)
    
    Ml <- as.numeric(rep(label,c) == rep(t(label_unique),each=n))
    Ml <- matrix(Ml,ncol=c)
    Mr <- as.numeric(rep(result,c) == rep(t(label_unique),each=n))
    Mr <- matrix(Mr,ncol=c)
    Pl<-colSums(Ml)/n
    Pr<-colSums(Mr)/n
    
    eps<-1e-10
    Hl <- -sum( Pl * log2( Pl + eps ) )
    Hr <- -sum( Pr * log2( Pr + eps ) )
    M <- t(Ml)%*%Mr/n
    Hlr <- -sum( M * log2(M + eps ) )
    
    MI <- Hl + Hr - Hlr
    v <- sqrt((MI/Hl)*(MI/Hr))
    return(v)
}



#' @title integrate evalution index
#' @description This fuction is used to integrate evaluation index.
#' @param c1 first label
#' @param c2 second label
#' @return
#' \item{res}{list of evaluation results}
#' @keywords evaluation
#' @examples 
#' \dontrun{
#' res<-evalt(c1,c2)
#' }
#' @export

evalt<-function(c1,c2){
    res1<-rand.idx(c1,c2)
    res2<-nmi(c1,c2)
    res3<-compute.f(c1,c2)
    res<-list(AR=res1$AR,RI=res1$RI,NMI=res2,f=res3$f,p=res3$p,r=res3$r)
    return(res)
}



#' @title compute positive part of matrix
#' @description This fuction is used to compute positive part of matrix.
#' @param A input matrix
#' @return
#' \item{Ap}{positive part of input matrix}
#' @keywords positive
#' @examples 
#' \dontrun{
#' Ap<-pos(A)
#' }
#' @export

pos<-function(A){
    Ap <- (A>=0)*A
    return(Ap)
}



#' @title compute negative part of matrix
#' @description This fuction is used to compute negative part of matrix.
#' @param A input matrix
#' @return
#' \item{Am}{negative part of input matrix}
#' @keywords negative
#' @examples 
#' \dontrun{
#' Am<-neg(A)
#' }
#' @export

neg<-function(A){
    Am <- (A<0)*(-A)
    return(Am)
}



#' @title run NNDSVD for initialization
#' @description This fuction is used to run NNDSVD for initialization.
#' @param A input matrix
#' @param k rank
#' @param flag mode for NNDSVD
#' @return
#' \item{W}{initial W for following NMF}
#' \item{H}{initial H for following NMF}
#' @keywords NNDSVD
#' @examples 
#' \dontrun{
#' ini<-NNDSVD(data,30,0)
#' }
#' @export

NNDSVD<-function(A,k,flag){

    W<-matrix(0,nrow=nrow(A),ncol=k)
    H<-matrix(0,nrow=k,ncol=ncol(A))
    res<-svd(A,nu=k,nv=k)
    U<-res$u
    V<-res$v
    S<-res$d
    W[,1]<-sqrt(S[1]) * abs(U[,1])      
    H[1,]<-sqrt(S[1]) * abs(t(V[,1])) 
    for(i in 2:k){
        uu<-U[,i] 
        vv<-V[,i]
        uup<-pos(uu)
        uun<-neg(uu) 
        vvp<-pos(vv)
        vvn<-neg(vv)
        n_uup<-norm(uup,'2')
        n_vvp<-norm(vvp,'2') 
        n_uun<-norm(uun,'2') 
        n_vvn<-norm(vvn,'2') 
        termp <- n_uup*n_vvp
        termn <- n_uun*n_vvn
        if (termp >= termn){
            W[,i] <- sqrt(S[i]*termp)*uup/n_uup
            H[i,] <- sqrt(S[i]*termp)*t(vvp)/n_vvp
        }
        else{
            W[,i] <- sqrt(S[i]*termn)*uun/n_uun 
            H[i,] <- sqrt(S[i]*termn)*t(vvn)/n_vvn
        }
            
    }
    W[which(W<0.0000000001)]<-0.1
    H[which(H<0.0000000001)]<-0.1
    if (flag==1){
        ind1<-which(W==0) 
        ind2<-which(H==0) 
        average<-mean(A) 
        W[ind1]<-average     
        H[ind2]<-average    
    }
    else if(flag==2){
        ind1<-which(W==0) 
        ind2<-which(H==0) 
        n1<-length(ind1)
        n2<-length(ind2)
        
        average<-mean(A)
        W[ind1]<-(average*stats::runif(n1,1)/100)  
        H[ind2]<-(average*stats::runif(n2,1)/100)     
    }
    return(list(W=W,H=H))
}



#' @title select rank
#' @description This fuction is used to select rank.
#' @param data data matrix
#' @return
#' \item{r}{selected rank}
#' @keywords rank
#' @examples 
#' \dontrun{
#' k<-select.r(data)
#' }
#' @export

select.r<-function(data){
    s<-svd(data)
    diags<-s$d
    diags<-diags-min(diags)
    ratiodiff<-rep(0,length(diags)-1)
    for(i in 1:(length(diags)-1)){
        ratiodiff[i]<-diags[i]/diags[i+1]
    }
    r<-which.min(ratiodiff<1.05)
    return(r)
}



#' @title compute soft threshold
#' @description This fuction is used to compute soft threshold.
#' @param A data matrix
#' @param gamma tuning parameters
#' @return
#' \item{S}{matrix output}
#' @keywords soft threshold
#' @examples 
#' \dontrun{
#' A<-soft.threshold(data,1.5)
#' }
#' @export

soft.threshold<-function(A,gamma){
    win1<-(A>0)
    win2<-(A<0)
    win3<-((gamma-2*abs(A))<0)
    S<-(A-0.5*gamma)*win1*win3+(A+0.5*gamma)*win2*win3
    return(S)
}



#' @title CRNMF for imputation
#' @description This fuction is used to run CRNMF for imputation.
#' @param data data matrix
#' @param maxitr max iteration number
#' @param lambda tuning parameter
#' @param r rank for NMF
#' @return
#' \item{W}{output W of CRNMF}
#' \item{H}{output H of CRNMF}
#' \item{S}{output imputation matrix S of CRNMF}
#' @keywords CRNMF, imputation
#' @examples 
#' \dontrun{
#' mat<-CRNMF(data,30,1.5,4)
#' }
#' @export

CRNMF<-function(data,maxitr=30,lambda=1.5,r){
    data<-data[which(rowSums(data)>0),]
    window<-(data==0) 
    S<-0
    res<-NNDSVD(data,r,0)
    W<-res$W
    H<-res$H
    cost_old<-2
    cost<-sum(sum(((data+S)-W%*%H)*((data+S)-W%*%H)))
    omega<-diag(colSums(data)/stats::median(colSums(data)))
    
    for(iter in 1:maxitr){
       k<-0
       while(abs(cost_old-cost)/(nrow(data)*ncol(data))>10^(-4)){
           k<-k+1
           W <- W*(((data+S)%*%t(H))/((W%*%H)%*%t(H))+1e-10)    
           W <- W/rep(sqrt(colSums(W^2)),each=nrow(W))
           H <- H*((t(W)%*%(data+S))/((t(W)%*%W)%*%H)+1e-10)
           cost_old<-cost
           cost<-sum(sum(((data+S)-W%*%H)*((data+S)-W%*%H)))
       }
       temp <- -data+W%*%H
       S <-soft.threshold(temp,rep(lambda*diag(omega),each=nrow(data)))
       S<-S*window
       cost_old<-cost
       cost<-sum(sum(((data+S)-W%*%H)*((data+S)-W%*%H)))
    }
    return(list(W=W,H=H,S=S))
}


