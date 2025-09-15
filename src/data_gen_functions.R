#'  Generate random params
#'
#' @description Add detail
#' @param Nrespondents Nrespondents
#' @param Qs Qs
#' @param Xs Xs
#' @importFrom stats rnorm
generate_params=function(Nrespondents,Qs,Xs){
  Nq_list=map(Qs,~dim(.)[1])
  Nq_total=sum(unlist(Nq_list))
  Nskill=dim(Qs[[1]])[2]
  Ntime=length(Qs)
  Nprofile=2^Nskill
  Ntransition=2^Ntime

  delta_cat=do.call(rbind,map(Qs,Q_to_delta))

  beta_mat=matrix(rnorm(Nprofile*Nq_total,sd=2),Nq_total,Nprofile)*
    delta_cat
  gamma_list=map(Xs,~map(.,function(x) rnorm(dim(x)[2])))
  return(list(beta_mat=beta_mat,gamma_list=gamma_list))
}

#'  Sample alpha from gamma
#'
#' @description Add detail
#' @param Nrespondents Nrespondents
#' @param Ntime Ntime
#' @param gamma_list gamma_list
#' @param Xs Xs
sample_alpha_from_gamma=function(Nrespondents,Ntime,gamma_list,Xs){
  trans_probs=gamma_to_transprobs(gamma_list,Xs)
  Nskill=length(Xs)
  trans_mat=matrix(NA,Nrespondents,Nskill)
  for(i in 1:Nrespondents){
    for(j in 1:Nskill){
      trans_mat[i,j]=sample(1:(2^Ntime),1,prob=trans_probs[[j]][i,])
    }
  }
  transitions=map_dfr(1:(2^Ntime)-1,
                      ~data.frame(t(rev(as.integer(intToBits(.))[1:Ntime]))))
  transitions_to_alpha(trans_mat,transitions)
}

#'  Generate data from betas and alphas
#'
#' @description Add detail
#' @param Nrespondents Nrespondents
#' @param Nq_list Nq_list
#' @param beta_mat beta_mat
#' @param alpha_mat alpha_mat
generate_data=function(Nrespondents,Nq_list,beta_mat,alpha_mat){
  Nprofile=dim(beta_mat)[2]
  profiles=map_dfr(1:Nprofile-1,
            ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
  A=model.matrix(data=profiles,as.formula(paste0('~', paste0(names(profiles),collapse='*'))))

  ltheta=beta_mat%*%t(A)
  theta=logistic(ltheta)
  Nq_total=sum(unlist(Nq_list))
  qt_map=unlist(map(1:length(Nq_list),~rep(.,Nq_list[[.]])))

  Ys=map(Nq_list,~matrix(NA,Nrespondents,.))
  jind_map=unlist(map(Nq_list,~1:.))
  for(i in 1:Nrespondents){
    for(j in 1:Nq_total){
      t=qt_map[j]
      p=theta[j,alpha_mat[i,t]]
      Ys[[t]][i,jind_map[j]]=sample(c(0,1),1,prob=c(1-p,p))
    }
  }
  Ys
}
