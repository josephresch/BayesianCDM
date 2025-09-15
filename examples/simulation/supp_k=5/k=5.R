{
  #gen data
  rm(list = ls())
  set.seed(300)
  setwd("~/Desktop/CDM/")
  devtools::load_all()
  
  require(mvtnorm)
  require(truncnorm)
  require(tidyverse)
  require(torch)
  library(doRNG)
  library(doParallel)
  
  
  
  M = 3000
  Nrespondents = 1000
  Nskill=k=5
  Ntime=2
  Nprofile=2^Nskill
  
  
  profiles=map_dfr(1:Nprofile-1,
                   ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
  profile_list = as.list(data.frame(t(profiles)))
  # Alternative for "profiles"
  # profile_list=list()
  # for(ii in 1:(2^Nskill)){
  #   tmp=rep(NA,Nskill)
  #   for(jj in 1:Nskill){
  #     tmp[Nskill-jj+1]=((ii-1)%/%(2^(jj-1)))%%2
  #   }
  #   profile_list[[ii]]=rev(tmp)
  # }
  
  
  
  # Q Matrix
  # Q=do.call(rbind,lapply(profile_list[-1],function(x)
  #   do.call(rbind,lapply(1:ceiling(16/(Nprofile-1)), function(i) x))))
  
  # No Interactions Q
  q_rep = 6
  Q = matrix(c(rep(c(0,0,0,0,1),q_rep),
               rep(c(0,0,0,1,0),q_rep),
               rep(c(0,0,1,0,0),q_rep), 
               rep(c(0,1,0,0,0),q_rep),
               rep(c(1,0,0,0,0),q_rep)),
             q_rep*Nskill, Nskill, byrow=T)
  
  # Random Q
  # Q=map(Nq_list,~profiles[sample(2^Nskill,Nq,replace=T),]) 
  
  # Manual Q
  # q_rep = 3
  # Q = matrix(c(rep(c(0,0,0),q_rep),
  #              rep(c(0,0,1),q_rep),
  #              rep(c(0,1,0),q_rep),
  #              rep(c(1,0,0),q_rep),
  #              rep(c(1,1,0),q_rep),
  #              rep(c(1,0,1),q_rep),
  #              rep(c(0,1,1),q_rep),
  #              rep(c(1,1,1),q_rep)), q_rep*8, 3, byrow=T)
  Qs      = list(Q,Q)
  Nq_list = map(Qs,~dim(.)[1])
  Nq      = dim(Q)[1]
  Nq_total=sum(unlist(Nq_list))
  
  
  
  
  # Standard one-intervention covariate
  # Xbase=cbind(rep(1,Nrespondents),c(rep(0,Nrespondents/2),rep(1,Nrespondents/2)))
  # Xs=build_Xlist_01(Xbase,Nskill,Ntime)
  
  # Complex covariates
  Intervention = c(rep(0, Nrespondents / 2), rep(1, Nrespondents / 2))
  Gender = sample(0:1, Nrespondents, replace = TRUE)
  Grade = sample(0:1, Nrespondents, replace = TRUE)
  ESL = sample(0:1, Nrespondents, replace = TRUE)
  SES = sample(0:1, Nrespondents, replace = TRUE)
  Xbase <- cbind(
    rep(1, Nrespondents),  # intercept
    Intervention,
    Gender,
    Grade,
    ESL,
    SES)
  Xs = build_Xlist_01(Xbase, Nskill, Ntime)
  
  
  delta=Q_to_delta(Q)
  true_params=list()
  true_params$beta_mat=
    cbind(rtruncnorm(Nq,mean=-1.5,sd=1,b=0),rtruncnorm(Nq,mean=3,sd=1,a=0),
          rtruncnorm(Nq,mean=2,sd=1,a=0),rtruncnorm(Nq,mean=2,sd=1,a=0),
          rtruncnorm(Nq,mean=2,sd=1,a=0), rtruncnorm(Nq,mean=2,sd=1,a=0),# intercept + base
          matrix(0,Nq,26))*   # interactions
    delta
  
  set.seed(1)
  gamma_means = c(0.5, 1, 2, -0.5, 1.5, 2.5)
  gamma_sds   = rep(0.1,6)
  true_params$gamma_list = purrr::map(Xs, function(skill_block) {
    purrr::map(skill_block, function(x) {
      p = ncol(x)
      rnorm(p, mean = gamma_means[1:p], sd = gamma_sds[1:p])
    })
  })
  
  exp(unlist(true_params$gamma_list))
  
  gamma_to_transprobs(true_params$gamma_list, Xs) 
  map(gamma_to_transprobs(true_params$gamma_list, Xs) , ~apply(.,2,mean))
}
plot(true_params$beta_mat[delta==1])


# Init params for covariates
# init_gamma = purrr::map(Xs, function(skill_block) {
#   purrr::map(skill_block, function(x) {rep(1, ncol(x)) })})
# init_beta = delta*2
# init_beta[,1] = -2
# init_params = list(beta_mat = init_beta, gamma_list = init_gamma)

init_params = NULL
# init_params = list(beta_mat = true_params$beta_mat, gamma_list = true_params$gamma_list)







# a=c(.7,.05,.05,.2)
# b=c(.2,.55,.05,.2)
# g1=log(a[2:4]/(1-sum(a[2:4])))
# g2=log(b[2:4]/(1-sum(b[2:4])))-g1
# gamma_list_true=map(1:Nskill,~
#   rbind(g1,g2),
# )
# Xtmp=build_Xlist_01(rbind(c(1,0),c(1,1)),Nskill,Ntime)
# gamma_list_true=
#   list(list(c(-2,3),c(-1.5),c(-1.5)),
#        list(c(-2,3),c(-1.5),c(-1.5)),
#        list(c(-2,3),c(-1.5),c(-1.5)))
# gamma_list=gamma_list_true
# true_params$gamma_list=gamma_list
# 
# gamma_to_transprobs(gamma_list,Xtmp)

# sim_results = list()

set.seed(1)
registerDoParallel(10)
Nsims=20


sim_results = foreach(sim = 1:Nsims) %dopar% {
  
  group_is=Xs[[1]][[1]][,2]+1
  alpha_sampdf=sample_alpha_from_gamma(Nrespondents,Ntime,true_params$gamma_list,Xs)%>%
    data.frame()%>%
    set_names(paste0('a',1:Ntime))%>%
    mutate(group=group_is)%>%
    relocate(group,.before=a1)%>%
    group_by(group)%>%
    arrange(group,a1,a2)
  true_alpha=alpha_sampdf%>%
    as.matrix()%>%
    .[,-1]
  # for additional covariates:
  true_alpha = sample_alpha_from_gamma(Nrespondents,Ntime,true_params$gamma_list,Xs)
  
  
  # simulate data
  profiles=map_dfr(1:Nprofile-1,
                   ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
  A=model.matrix(data=profiles,as.formula(paste0('~', paste0(names(profiles),collapse='*'))))
  
  ltheta=true_params$beta_mat%*%t(A)
  theta=logistic(ltheta)
  qt_map=unlist(map(1:length(Nq_list),~rep(.,Nq_list[[.]])))
  
  Ys=map(Nq_list,~matrix(NA,Nrespondents,.))
  for (t in 1:Ntime){
    for(i in 1:Nrespondents){
      for(j in 1:Nq){
        p=theta[j,true_alpha[i,t]]
        Ys[[t]][i,j]=sample(c(0,1),1,prob=c(1-p,p))
      }
    }
  }
  
  priors=list(beta_prior=2.5,gamma_prior=1)
  
  runtime = system.time({sampler_out= fit_longitudinal_cdm_full(Ys,Xs,Qs,M,
                                          initparams = if(is.null(init_params)) NULL else init_params,
                                          priors=priors,
                                          fixed_beta = TRUE)})
  # rm burn-in
  sampler_out$samples_burned_in = sampler_out$samples[500:M,]
  
  list(sim = sim, true_alpha = true_alpha, Ys = Ys, 
       sampler_out = sampler_out)
}


#################################################################################################################################################################

{
  beta_list = lapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                       dplyr::select(contains('beta')))
  full_beta_mat = do.call(rbind, beta_list)
  # boxplot(full_beta_mat[,60:78])
  true_alphas = lapply(1:Nsims, function(x) sim_results[[x]]$true_alpha)
  Ys = lapply(1:Nsims, function(x) sim_results[[x]]$Ys)
  # Beta cleaning
  beta_mat=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                      dplyr::select(contains('beta'))%>%
                      apply(2,mean)))
  beta_sd=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                     dplyr::select(contains('beta'))%>%
                     apply(2,sd)))
  # Gamma cleaning
  gamma_post=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                        dplyr::select(contains('gamma'))%>%
                        apply(2,mean)))
  gamma_sd=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                      dplyr::select(contains('gamma'))%>%
                      apply(2,sd)))
  Nmain = sum(delta[,2:(Nskill+1)])
  Nint = sum(delta[,(Nskill+2):ncol(delta)])
}


plot(true_params$beta_mat[delta==1])
### Beta posteriors
ggplot(data.frame(beta.post=c(beta_mat), q=rep(1:sum(delta), each = Nsims), 
                  beta.true=rep(true_params$beta_mat[delta==1], each=Nsims),
                  term=factor(c(rep("Intercepts", Nq*Nsims), rep("Main Effects", Nmain*Nsims), 
                                rep("Interactions", Nint*Nsims)), 
                              levels=c("Intercepts","Main Effects","Interactions")))) + 
  scale_color_manual(values = c("#56B4E9", "#E69F00", "seagreen3")) +
  geom_point(aes(x=q ,y=beta.true, col=term), size=.7)+  
  geom_boxplot(aes(y=beta.post, x=q, group=q, col = term), outlier.size=0.05,
               fill = "white", position="identity", alpha=.2, width=1, size=0.4)+
  theme(legend.position = c(.85,.15),
        legend.background = element_rect(fill="linen",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="white"),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17))+
  xlab("") +
  ylab(expression(Beta ~"Value")) +
  labs(color='Posterior Distributions \n+ True Values:') 



### Gamma posteriors
gdim = 8
plotdf=data.frame(gamma_post=c(gamma_post),
                  gamma_true=rep(unlist(true_params$gamma_list),each=Nsims),
                  i=1:length(c(gamma_post)))%>%
  mutate(skill=c(rep("Attribute 1",gdim*Nsims),rep("Attribute 2",gdim*Nsims),
                 rep("Attribute 3",gdim*Nsims),rep("Attribute 4",gdim*Nsims), 
                 rep("Attribute 5",gdim*Nsims)),
         i=factor(rep(c(1:gdim,1:gdim,1:gdim,1:gdim,1:gdim),each=Nsims)), 
         q=rep(1:(k*gdim),each=Nsims))
# ggplot
a=ggplot(plotdf)+
  geom_boxplot(aes(y=gamma_post, x=i, group=q, col='Posterior'), outlier.size=0,
               fill = "white", position="identity", alpha=.5)+  
  geom_point(aes(x=i,y=gamma_true,col='True'))+
  scale_x_discrete(guide=guide_axis(angle = 90),
                   labels=c(expression(atop(NA, atop(0%->%1,"Intercept"))),
                            expression(atop(NA, atop(0%->%1,"Covariate 1"))),
                            expression(atop(NA, atop(0%->%1,"Covariate 2"))),
                            expression(atop(NA, atop(0%->%1,"Covariate 3"))),
                            expression(atop(NA, atop(0%->%1,"Covariate 4"))),
                            expression(atop(NA, atop(0%->%1,"Covariate 5"))),
                            expression(atop(NA, atop(1%->%0,"Intercept"))),
                            expression(atop(NA, atop(1%->%1,"Intercept"))))) +
  facet_grid(~skill) +
  theme(legend.position = c(0.77, 0.12),
        legend.background = element_rect(fill="linen",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="white"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(size = 17))+
  labs(color=NULL) +
  xlab("") +
  ylab(expression(Gamma ~ "Value")) 
ggsave(a,file='~/Desktop/CDM/plots/k5_sim_gamma.pdf',height=6,width=8)




# Beta credible interval coverage
beta_upper = beta_mat + 2*beta_sd
beta_lower = beta_mat - 2*beta_sd
true_beta = true_params$beta_mat[delta==1]
beta_credible_post = matrix(0,Nsims,length(true_beta))
for(i in 1:length(true_beta))
{
  beta_credible_post[,i] = (beta_lower[,i] < true_beta[i]) & (true_beta[i] < beta_upper[,i])
}
data.frame(prob = round(apply(beta_credible_post,2,mean),3),
           i = 1:length(true_beta),
           term = factor(c(rep("Intercepts", Nq), rep("Main Effects", Nmain), 
                           rep("Interactions", Nint)), 
                         levels=c("Intercepts","Main Effects","Interactions"))) %>%
  ggplot(aes(y=prob, x=i, group = term, col=term))+
  geom_bar(stat="identity") +
  scale_color_manual(values = c("seagreen3", "#E69F00", "#56B4E9")) +
  theme(legend.position = c(0.8, 0.2),
        legend.background = element_rect(fill="linen",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="white"),
        text = element_text(size = 17),
        legend.text = element_text(size = 12))+
  xlab("") +
  ylab("Coverage Probability") +
  labs(color=expression(Beta ~"Coefficients")) 

# Gamma credible interval coverage
gamma_upper = gamma_post + 2*gamma_sd
gamma_lower = gamma_post - 2*gamma_sd
true_gamma = unlist(true_params$gamma_list)
gamma_credible_post = matrix(0,Nsims,length(true_gamma))
for(i in 1:length(true_gamma))
{
  gamma_credible_post[,i] = (gamma_lower[,i] < true_gamma[i]) & (true_gamma[i] < gamma_upper[,i])
}
data.frame(prob = round(apply(gamma_credible_post,2,mean),3),
           i = 1:length(true_gamma),
           term = factor(c(rep("Skill 1", length(true_gamma)/Nskill), 
                           rep("Skill 2", length(true_gamma)/Nskill), 
                           rep("Skill 3", length(true_gamma)/Nskill), 
                           rep("Skill 4", length(true_gamma)/Nskill), 
                           rep("Skill 5", length(true_gamma)/Nskill)), 
                         levels=c("Skill 1","Skill 2","Skill 3", 
                                  "Skill 4", "Skill 5"))) %>%
  ggplot(aes(y=prob, x=i, group = term, col=term))+
  geom_bar(stat="identity") +
  scale_color_manual(values = c("seagreen3", "#E69F00", "#56B4E9",
                                "white", "black")) +
  theme(legend.position = c(0.8, 0.2),
        legend.background = element_rect(fill="linen",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="white"),
        text = element_text(size = 17),
        legend.text = element_text(size = 12))+
  xlab("") +
  ylab("Coverage Probability") +
  labs(color=expression(Gamma ~"Coefficients"))



