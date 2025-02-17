The zip file contains the following files: 

1)mediation.r with the function mediation(); using the command source("mediation.r") to install. Once  
         mediation(ppt,data,int="TRUE",PX=NULL,tol=0.005,setting,level=lv): main function to conduct the mediation estimation; 
                 ppt: a scale or a vector, input of the user-interested time of the time-related functions: Omega_DE, Omega_IE, and log pi_Omegq(page 4 of main text)
                data: data.frame format:
                            1st column: event times
		                2nd column: indicator of the event, 0: censored, 1: the event occurs
		                3rd column: indicator of the selection, 0: unselected, 1: selected
		                4th column: binary exposure 0/1
		                5th column: value of the mediator
                        6th-... column: covariates
                int: with/without ("TRUE"/"FALSE") considering the interaction effect in the model     
                 PX: The vector of chosen columns, the columns are treated as the covariates to calculate the weight (see the model in simulation of the main text on page 10) 
                tol: the desired accuracy (convergence tolerance) with default is 0.005
            setting: setting the level of the covariate at thire mean/median/0/user_setting ("center"/"median"/"zero"/"level")
              level: lv (scale or vector),  the level of covariates by the user setting
Please see the explanation of the output in the file: mediation.r
 
1-1)probit.r with the functions predit.haz() and tran.npmle(); using the command source("probit.r") to install(Ver4).
        predit.haz(): sub-function(tool) in the mediation()  
        tran.npmle(): sub-function(tool) in the mediation()  
    
  
2)IV.r with the function IV(); using the command source("IV.r") to install. 
 IVfunction(ppt,data,m0=0,m1=1,PX=NULL,tol=0.005,SS=1,mu=0,bin=0.5,setting,level=lv): main function to conduct the estimation of IV analyses; 
                 ppt: a scale or a vector, input of the user-interested time of the time-related functions: Omega_DE, Omega_IE, and log pi_Omegq(page 4 of main text)
                data: data.frame format:
                            1st column: event times
		                2nd column: indicator of the event, 0: censored, 1: the event occurs
		                3rd column: indicator of the selection, 0: unselected, 1: selected
		                4th column: binary exposure 0/1
		                5th column: value of the mediator
                        6th-... column: covariates
             m0(m1): level of the mediator 
                 PX: The vector of chosen columns, the columns are treated as the covariates to calculate the weight (see the model in simulation of the main text on page 10) 
                tol: the desired accuracy (convergence tolerance) with default is 0.005
                 SS: the standard deviation of error term  e_Ti (page 6)
            setting: setting the level of the covariate at thire mean/median/0/user_setting ("center"/"median"/"zero"/"level")
              level: lv (scale or vector),  the level of covariates by the user setting
       Please see the explanation of the output in the file: IV.r
2-1)probit_W.r with the function tran.npmle(); using the command source("probit_W.r") to install(ver5). 
        tran.npmle(): sub-function(tool) in the IV 13)




Binding antibody N.RData: an R data containing the real data as a data.frame (consider the exposure as binding antibody N)
Binding antibody RBD.RData: an R data containing the real data as a data.frame (consider the exposure as binding antibody N)
Binding antibody Spike.RData: an R data containing the real data as a data.frame (consider the exposure as binding antibody N)
Pseudovirus 50.RData: an R data containing the real data as a data.frame (consider the exposure as pseudovirus 80)
Pseudovirus 80.RData: an R data containing the real data as a data.frame (consider the exposure as pseudovirus 80)
Wild Type Live Virus.RData: an R data containing the real data as a data.frame (consider the exposure as wild type live virus)



example_Real data.R: an example code of applying the proposed mediation tests to the REVEAL study    



                                        