
########################################
## GENERAL PACKAGE DOCUMENTATION PAGE ##
########################################


#' Estimating Disease Epidemiology
#'
#' DisModR is a version of DisMod II which estimates disease incidence, remission, case-fatality, mortality, and prevalence using as input any of the three expected outputs.
#'
#' @docType package
#' @name DisModR
NULL


########################################
## IMPORTING OTHER PACKAGES/FUNCTIONS ##
########################################


#' @import graphics grDevices stats neldermead
#' @importFrom graphics plot legend lines
#' @importFrom grDevices recordPlot
#' @importFrom stats smooth.spline predict
#' @importFrom neldermead fminsearch
NULL

###################
## MAIN FUNCTION ##
###################

#' @export
#' @title Estimating Disease Epidemiology through Analytical or Optimization Methods.
#' @description This function provides estimates and graphs for the epidemiology of a single disease.
#' @param Data Data object containing the input variables. Contains the columns for
#' AgeLB, AgeUB, Incidence, Remission, Case_Fatality, Mortality, Prevalence, and/or Total_Mortality_Group.
#' @param AgeLB Lower bound for age groups
#' @param AgeUB Upper bound for age groups
#' @param Incidence Incidence hazard of selected disease
#' @param Remission Remission hazard of selected disease
#' @param Case_Fatality Case Fatality hazard of selected disease
#' @param Mortality Disease-specific mortality hazard of the population
#' @param Prevalence Prevalence proportion of the selected disease in the population
#' @param RR_Mortality Relative risk mortality of the selected disease
#' @param Total_Mortality Total mortality rates of a population. Use this if data for total mortality
#' is in one-year age intervals. In this case, the input must come from a separate object other than the Data list.
#' @param Total_Mortality_Group Total mortality rates of a population. Use this if data for total mortality
#' is in group intervals like the other inputs. In this case, the input must come from the Data list.
#' @param Cubic_spline Use cubic splines to smooth the grouped data. Default is TRUE.
#' @param birth_prevalence Birth prevalence. Prevalence at age 0. Default is 0.
#' @param spar If Cubic_spline is TRUE, can be set manually by the user. Smoothing parameter,
#' typically (but not necessarily) in (0,1]. Default is NULL. See \emph{smooth.spline} for more details.
#' @param cv if Cubic_spline is TRUE, can be set manually by the user. Ordinary leave-one-out (TRUE)
#' or ‘generalized’ cross-validation (GCV) when FALSE; is used for smoothing parameter computation only when
#' spar is not specified. Default is TRUE. See \emph{smooth.spline} for more details.
#' @param Tol An optimization option which drives the behavior of \emph{fminsearch}.  It
#' is the absolute tolerance on simplex size. The default value is 1.0E-9.
#' @param MaxIt An optimization option which drives the behavior of \emph{fminsearch}. It
#' is the maximum number of iterations. The default is 1000.
#' @return A nested list containing the list for the inputs, cubic-spline adjusted inputs, final estimates in yearly intervals,
#' final estimates in their group intervals, as well as their corresponding graph functions.
#' Use $Input to call the Input table (converted to yearly intervals), $Plot_input to cal
#' the plot for the input values, $Cubic_spline to call for the input data smoothed using
#' cubic splines, $Plot_cubic_spline to call for the plot of the Cubic spline smoothed data,
#' $Final_estimates to call for the final estimates in yearly form, $Plot_Final_estimates
#' to call for the plot of the Final estimates in yearly form, $Final_Table to call for
#' the final estimates in group intervals (as in the established groups in the inputs), and
#' $Plot_Final_groupestimates to call for the plot of the final estimates in group intervals.
#' Other outputs such as RR_Mortality can be obtained by calling the list name $RR_Mortality
#' within the $Final_estimates list.
#' @author Jade Vincent Membrebe, Ziv Shkedy, Brecht Devleesschauwer, Ewoud De Troyer
#' @source This function is based on the program DisMod II.
#' @references Barendregt, J. J., Van Oortmarssen, G. J., Vos, T., & Murray, C. J. (2003).
#' A generic model for the assessment of disease epidemiology: the computational basis
#' of DisMod II. \emph{Population health metrics, 1(1)}, 4.


DisMod <- function (Data, AgeLB, AgeUB, Incidence=NULL, Remission=NULL,
                    Case_Fatality=NULL,Mortality=NULL, Prevalence=NULL, RR_Mortality=NULL,
                    Total_Mortality=NULL, Total_Mortality_Group=NULL, birth_prevalence=0,
                    Cubic_spline=TRUE, spar=NULL, cv=TRUE, Tol=1.0E-9, MaxIt=1000){

  ###########################################################################
  ##########                     Input Function                    ##########
  ###########################################################################

  Epi_data <- function (AgeLB, AgeUB, Incidence, Remission,
                        Case_Fatality, Mortality, Prevalence, RR_Mortality,
                        Total_Mortality, Total_Mortality_Group){

    Age1=NULL
    Age2=NULL
    Incidence1=NULL
    Incidence2=NULL
    Remission1=NULL
    Remission2=NULL
    Case_Fatality1=NULL
    Case_Fatality2=NULL
    Mortality1=NULL
    Mortality2=NULL
    Prevalence1=NULL
    Prevalence2=NULL
    RR_Mortality1=NULL
    RR_Mortality2=NULL
    Total_Mortality1=NULL
    Total_Mortality2=NULL

    for (i in 1:length(AgeUB)){
      for (j in 1:AgeUB[i]){
        if(j>=AgeLB[i]){
          Age1[j]=j
          Incidence1[j]=Incidence[i]
          Remission1[j]=Remission[i]
          Case_Fatality1[j]=Case_Fatality[i]
          Mortality1[j]=Mortality[i]
          Prevalence1[j]=Prevalence[i]
          RR_Mortality1[j]=RR_Mortality[i]
          Total_Mortality1[j]=Total_Mortality_Group[i]
        }}}

    for (ag in 1:length(Age1)){
      Age2[ag+1]=Age1[ag]
      Incidence2[ag+1]=Incidence1[ag]
      Remission2[ag+1]=Remission1[ag]
      Case_Fatality2[ag+1]=Case_Fatality1[ag]
      Mortality2[ag+1]=Mortality1[ag]
      Prevalence2[ag+1]=Prevalence1[ag]
      RR_Mortality2[ag+1]=RR_Mortality1[ag]
      if (!is.null(Total_Mortality1)){
        Total_Mortality2[ag+1]=Total_Mortality1[ag]
      } else {Total_Mortality2[ag+1]=Total_Mortality[ag]}
    }

    Age2[1]=0
    if (!is.null(Incidence2)){Incidence2[1]=Incidence2[2]}
    if (!is.null(Remission2)){Remission2[1]=Remission2[2]}
    if (!is.null(Case_Fatality2)){Case_Fatality2[1]=Case_Fatality2[2]}
    if (!is.null(Mortality2)){Mortality2[1]=Mortality2[2]}
    if (!is.null(Prevalence2)){Prevalence2[1]=Prevalence2[2]}
    if (!is.null(RR_Mortality2)){RR_Mortality2[1]=RR_Mortality2[2]}
    if (!is.null(Total_Mortality2)){Total_Mortality2[1]=Total_Mortality2[2]}

    output=list(Age=Age2,Incidence=Incidence2,Remission=Remission2,
                Case_Fatality=Case_Fatality2,Mortality=Mortality2,
                Prevalence=Prevalence2, RR_Mortality=RR_Mortality2,
                Total_Mortality=Total_Mortality2)
    class(output) <- "Epi_data"
    return(output)
  }

  plot.Epi_data <- function (Epi_raw){
    par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
    plot(Epi_raw$Age,Epi_raw$Age,type="n",col=4,lwd=3,
         ylim=c(min(c(Epi_raw$Incidence,Epi_raw$Remission,Epi_raw$Mortality,
                      Epi_raw$Case_Fatality,Epi_raw$Prevalence)),
                max(c(Epi_raw$Incidence,Epi_raw$Remission,Epi_raw$Mortality,
                      Epi_raw$Case_Fatality,Epi_raw$Prevalence))),
         ylab=c("Hazard or Rate"),xlab="Age",main="Input Variables")
    count=0
    legend("topright", inset=c(-0.2,count), legend=c(""), title="Legend:", bty="n")
    if (length(Epi_raw$Age)==length(Epi_raw$Incidence)){
      lines(Epi_raw$Age,Epi_raw$Incidence,col="blue",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.4,count),legend=c("Incidence"),lwd=3, col="blue", bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Remission)){
      lines(Epi_raw$Age,Epi_raw$Remission,col="yellow",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.42,count),legend=c("Remission"),lwd=3,col="yellow",bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Mortality)){
      lines(Epi_raw$Age,Epi_raw$Mortality,col="red",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.385,count),legend=c("Mortality"),lwd=3,col="red",bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Case_Fatality)){
      lines(Epi_raw$Age,Epi_raw$Case_Fatality,col="purple",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.458,count),legend=c("Case Fatality"),lwd=3,col="purple",bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Prevalence)){
      lines(Epi_raw$Age,Epi_raw$Prevalence,col="green",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.43,count),legend=c("Prevalence"),lwd=3,col="green",bty="n")}
  }

  Input=Epi_data(AgeLB=Data$AgeLB, AgeUB=Data$AgeUB, Incidence=Data$Incidence,
                 Remission=Data$Remission, Case_Fatality=Data$Case_Fatality,
                 Mortality=Data$Mortality, Prevalence=Data$Prevalence, RR_Mortality=Data$RR_Mortality,
                 Total_Mortality=Total_Mortality, Total_Mortality_Group=Data$Total_Mortality_Group)

  plot(Input)
  plot_input=recordPlot()

  ###########################################################################
  ##########                     Cubic Splines                     ##########
  ###########################################################################

  if (Cubic_spline==TRUE){
    Epi_data2 <- function (Data,spar=spar,cv=cv) {

      input_inci=NULL
      input_remi=NULL
      input_cfat=NULL
      input_mort=NULL
      input_prev=NULL
      input_rrmort=NULL


      if (!is.null(Data$Incidence)){
        Incidence_spline=smooth.spline(x=Data$Age,y=Data$Incidence,spar=spar,cv=cv)
        input_inci1=predict(Incidence_spline,Data$Age)
        input_inci=input_inci1[[2]]
        for (ag in 1:(length(Data$Age)-1)){
          input_inci[1]=0
          if (input_inci[ag+1]<0) {
            input_inci[ag+1]=input_inci[ag]
          }}
      }

      if (!is.null(Data$Remission)){
        Remission_spline=smooth.spline(x=Data$Age,y=Data$Remission,spar=spar,cv=cv)
        input_remi1=predict(Remission_spline,Data$Age)
        input_remi=input_remi1[[2]]
        for (ag in 1:(length(Data$Age)-1)){
          input_remi[1]=0
          if (input_remi[ag+1]<0) {
            input_remi[ag+1]=input_remi[ag]
          }}
      }

      if (!is.null(Data$Case_Fatality)){
        Case_Fatality_spline=smooth.spline(x=Data$Age,y=Data$Case_Fatality,spar=spar,cv=cv)
        input_cfat1=predict(Case_Fatality_spline,Data$Age)
        input_cfat=input_cfat1[[2]]
        for (ag in 1:(length(Data$Age)-1)){
          input_cfat[1]=0
          if (input_cfat[ag+1]<0) {
            input_cfat[ag+1]=input_cfat[ag]
          }}
      }

      if (!is.null(Data$Mortality)){
        Mortality_spline=smooth.spline(x=Data$Age,y=Data$Mortality,spar=spar,cv=cv)
        input_mort1=predict(Mortality_spline,Data$Age)
        input_mort=input_mort1[[2]]
        for (ag in 1:(length(Data$Age)-1)){
          input_mort[1]=0
          if (input_mort[ag+1]<0) {
            input_mort[ag+1]=input_mort[ag]
          }}
      }

      if (!is.null(Data$Prevalence)){
        Prevalence_spline=smooth.spline(x=Data$Age,y=Data$Prevalence,spar=spar,cv=cv)
        input_prev1=predict(Prevalence_spline,Data$Age)
        input_prev=input_prev1[[2]]
        for (ag in 1:(length(Data$Age)-1)){
          input_prev[1]=0
          if (input_prev[ag+1]<0) {
            input_prev[ag+1]=input_prev[ag]
          }}
      }

      if (!is.null(Data$RR_Mortality)){
        RR_Mortality_spline=smooth.spline(x=Data$Age,y=Data$RR_Mortality,spar=spar,cv=cv)
        input_rrmort1=predict(RR_Mortality_spline,Data$Age)
        input_rrmort=input_rrmort1[[2]]
        for (ag in 1:(length(Data$Age)-1)){
          input_rrmort[1]=0
          if (input_rrmort[ag+1]<0) {
            input_rrmort[ag+1]=input_rrmort[ag]
          }}
      }

      output <- list(Age=Data$Age,input_inci=input_inci,input_remi=input_remi,
                     input_cfat=input_cfat, input_mort=input_mort, input_prev=input_prev,
                     input_rrmort=input_rrmort, totmor=Data$Total_Mortality)
      class(output) <- "Epi_data2"
      return(output)
    }

    plot.Epi_data2 <- function(Epi_raw2){
      Epi_raw2=lapply(Epi_raw2, "[", -1)
      par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
      plot(Epi_raw2$Age,Epi_raw2$Age,lwd=3, col=4,type="n",
           ylim=c(min(c(Epi_raw2$input_inci,Epi_raw2$input_remi,Epi_raw2$input_cfat,
                        Epi_raw2$input_mort,Epi_raw2$input_prev)),
                  max(c(Epi_raw2$input_inci,Epi_raw2$input_remi,Epi_raw2$input_cfat,
                        Epi_raw2$input_mort,Epi_raw2$input_prev))),
           ylab=c("Hazard or Rate"),xlab="Age",main="Input Variables Adjusted by Cubic Splines")
      count=0
      legend("topright", inset=c(-0.2,count), legend=c(""), title="Legend:", bty="n")
      if (!is.null(Epi_raw2$input_inci)){
        lines(Epi_raw2$input_inci,lwd=3,col="blue")
        count=count+0.1
        legend("topright",inset=c(-0.4,count),legend=c("Incidence"),lwd=3, col="blue", bty="n")}
      if (!is.null(Epi_raw2$input_remi)){
        lines(Epi_raw2$input_remi,lwd=3,col="yellow")
        count=count+0.1
        legend("topright",inset=c(-0.42,count),legend=c("Remission"),lwd=3,col="yellow",bty="n")}
      if (!is.null(Epi_raw2$input_mort)){
        lines(Epi_raw2$input_mort,lwd=3,col="red")
        count=count+0.1
        legend("topright",inset=c(-0.385,count),legend=c("Mortality"),lwd=3,col="red",bty="n")}
      if (!is.null(Epi_raw2$input_cfat)){
        lines(Epi_raw2$input_cfat,lwd=3,col="purple")
        count=count+0.1
        legend("topright",inset=c(-0.458,count),legend=c("Case Fatality"),lwd=3,col="purple",bty="n")}
      if (!is.null(Epi_raw2$input_prev)){
        lines(Epi_raw2$input_prev,lwd=3,col="green")
        count=count+0.1
        legend("topright",inset=c(-0.43,count),legend=c("Prevalence"),lwd=3,col="green",bty="n")}
    }

    Cubic_spline_data=Epi_data2(Data=Input,spar=spar,cv=cv)
    Compartmental_model=Cubic_spline_data
    plot(Cubic_spline_data)
    plot_cubic_spline=recordPlot()
  } else {
    Epi_data2_1 <- function (Data) {

      input_inci=NULL
      input_remi=NULL
      input_cfat=NULL
      input_mort=NULL
      input_prev=NULL
      input_rrmort=NULL


      if (!is.null(Data$Incidence)){
        input_inci=Data$Incidence
      }

      if (!is.null(Data$Remission)){
        input_remi=Data$Remission
      }

      if (!is.null(Data$Case_Fatality)){
        input_cfat=Data$Case_Fatality
      }

      if (!is.null(Data$Mortality)){
        input_mort=Data$Mortality
      }

      if (!is.null(Data$Prevalence)){
        input_prev=Data$Prevalence
      }

      if (!is.null(Data$RR_Mortality)){
        input_rrmort=Data$RR_Mortality
      }

      output <- list(Age=Data$Age,input_inci=input_inci,input_remi=input_remi,
                     input_cfat=input_cfat, input_mort=input_mort, input_prev=input_prev,
                     input_rrmort=input_rrmort, totmor=Data$Total_Mortality)
      return(output)
    }
    Cubic_spline_data=NULL
    plot_cubic_spline=NULL
    Compartmental_model=Epi_data2_1(Data=Input)
  }

  if(!is.null(Compartmental_model$input_inci)&!is.null(Compartmental_model$input_remi)&
     !is.null(Compartmental_model$input_cfat)){

    ###########################################################################
    ##########  I Accidentally Vaporize My Maths Teacher             ##########
    ###########################################################################

    IRC <- function (input_inci, input_remi,input_cfat, totmor, Age, bprev=0){

      ###########################################################################
      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_prev=NULL
        input_dura=NULL
        input_mort=NULL
        input_rrmort=NULL
        incipoprate=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }
      ###########################################################################

      for (ag in 1:(length(Age)-1)){

        ###########################################################################
        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            input_mort[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {input_mort[ag+1]=input_mort[ag]}

          input_prev[ag+1]=0.5*(state_withcondition[ag]
                                +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }
        ###########################################################################

        input_mort[1]=0
        input_prev[1]=0
      }	#For-loop end#

      output=(list(Incidence=input_inci,Remission=input_remi,
                   Case_Fatality=input_cfat,Prevalence_proportion=input_prev,
                   Mortality_rate=input_mort,RR_Mortality=input_rrmort))
      class(output)="Epi_final"
      return(output)
    }


    ###########################################################################
    ###########################################################################
    ###########################################################################
    Final_estimates=IRC(input_inci=Compartmental_model$input_inci,
                        input_remi=Compartmental_model$input_remi,
                        input_cfat=Compartmental_model$input_cfat,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }


  if(!is.null(Compartmental_model$input_inci)&!is.null(Compartmental_model$input_remi)&
     !is.null(Compartmental_model$input_prev)){

    ###########################################################################
    ##########  Three Old Ladies Knit the Socks of Death             ##########
    ###########################################################################

    IRP <- function (input_inci, input_remi, input_prev, totmor, Age, bprev=0,
                     fminopt = optimset(TolX=Tol,MaxIter=MaxIt)){

      ###########################################################################
      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_mort=NULL
        input_dura=NULL
        input_rrmort=NULL
        incipoprate=NULL
        input_cfat=NULL
        input_cfat2=NULL
        input_prev2=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }
      ###########################################################################

      for (ag in 1:(length(Age)-1)){

        ###########################################################################
        amoeba <- function(x){

          ### intermediate variables
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+x[1]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*x[1]+input_remi[ag]**2+2.0*x[1]*input_remi[ag]+x[1]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (x[1]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((x[1]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(x[1]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          input_mort[ag+1]=(state_dead[ag+1]-state_dead[ag])/py

          (input_prev[ag+1]-0.5*(state_withcondition[ag]
                                 +state_withcondition[ag+1])/py)**2

        }


        neldermead=fminsearch(fun=amoeba, x0=neldtmp2[ag], fminopt)
        neldtmp1=neldermead$simplexopt
        neldtmp2[ag+1]=neldtmp1$x[1,]

        input_cfat[1]=0

        if (neldtmp2[ag+1]>=0 & abs(neldtmp2[ag+1]-neldtmp2[ag])<0.2)   {  #
          input_cfat[ag+1]=neldtmp2[ag+1]
        } else {input_cfat[ag+1]=input_cfat[ag]}
        ###########################################################################


        ###########################################################################
        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            intermediate_mortotaal[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {intermediate_mortotaal[ag+1]=intermediate_mortotaal[ag]}

          input_prev2[ag+1]=0.5*(state_withcondition[ag]
                                 +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }
        ###########################################################################

        intermediate_mortotaal[1]=0
        input_prev2[1]=0
      }	#For-loop end#
      cfat_spline=smooth.spline(x=Age,y=input_cfat,spar=0.4)
      input_cfat1=predict(cfat_spline,Age)
      input_cfat_downhill=input_cfat1[[2]]
      for (ag in 1:(length(Age)-1)){
        input_cfat_downhill[1]=0
        if (input_cfat_downhill[ag+1]<0) {
          input_cfat_downhill[ag+1]=input_cfat_downhill[ag]
        }}
      output=(list(Incidence=input_inci,Remission=input_remi,
                   Case_fatality_downhill=input_cfat_downhill,Prevalence_proportion=input_prev,
                   Mortality_rate=intermediate_mortotaal,RR_Mortality=input_rrmort))
      class(output)="Epi_final"
      return(output)
    }


    ###########################################################################
    ###########################################################################
    ###########################################################################
    Final_estimates=IRP(input_inci=Compartmental_model$input_inci,
                        input_remi=Compartmental_model$input_remi,
                        input_prev=Compartmental_model$input_prev,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }


  if(!is.null(Compartmental_model$input_inci)&!is.null(Compartmental_model$input_remi)&
     !is.null(Compartmental_model$input_mort)){

    ###########################################################################
    ##########  Grover Unexpectedly Loses His Trouser                ##########
    ###########################################################################

    IRM <- function (input_inci, input_remi, input_mort, totmor, Age, bprev=0,
                     fminopt = optimset(TolX=Tol, MaxIter=MaxIt)){

      ###########################################################################
      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_prev=NULL
        input_dura=NULL
        input_rrmort=NULL
        incipoprate=NULL
        input_cfat=NULL
        input_cfat2=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }
      ###########################################################################

      for (ag in 1:(length(Age)-1)){

        ###########################################################################
        amoeba <- function(x){

          ### intermediate variables
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+x[1]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*x[1]+input_remi[ag]**2+2.0*x[1]*input_remi[ag]+x[1]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (x[1]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((x[1]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(x[1]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          (input_mort[ag+1]-((state_dead[ag+1]-state_dead[ag])/py))**2
        }



        neldermead=fminsearch(fun=amoeba, x0=neldtmp2[ag], fminopt)
        neldtmp1=neldermead$simplexopt
        neldtmp2[ag+1]=neldtmp1$x[1,]

        input_cfat[1]=0

        if (neldtmp2[ag+1]>=0 & abs(neldtmp2[ag+1]-neldtmp2[ag])<0.2)   {  #
          input_cfat[ag+1]=neldtmp2[ag+1]
        } else {input_cfat[ag+1]=input_cfat[ag]}
        ###########################################################################


        ###########################################################################
        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            intermediate_mortotaal[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {intermediate_mortotaal[ag+1]=intermediate_mortotaal[ag]}

          input_prev[ag+1]=0.5*(state_withcondition[ag]
                                +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }
        ###########################################################################

        intermediate_mortotaal[1]=0
        input_prev[1]=0
      }	#For-loop end#
      cfat_spline=smooth.spline(x=Age,y=input_cfat,spar=0.4)
      input_cfat1=predict(cfat_spline,Age)
      input_cfat_downhill=input_cfat1[[2]]
      for (ag in 1:(length(Age)-1)){
        input_cfat_downhill[1]=0
        if (input_cfat_downhill[ag+1]<0) {
          input_cfat_downhill[ag+1]=input_cfat_downhill[ag]
        }}
      output=(list(Incidence=input_inci,Remission=input_remi,
                   Case_fatality_downhill=input_cfat_downhill,Prevalence_proportion=input_prev,
                   Mortality_rate=intermediate_mortotaal,RR_Mortality=input_rrmort))
      class(output)="Epi_final"
      return(output)
    }


    ###########################################################################
    ###########################################################################
    ###########################################################################
    Final_estimates=IRM(input_inci=Compartmental_model$input_inci,
                        input_remi=Compartmental_model$input_remi,
                        input_mort=Compartmental_model$input_mort,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }


  if(!is.null(Compartmental_model$input_inci)&!is.null(Compartmental_model$input_cfat)&
     !is.null(Compartmental_model$input_prev)){

    ###########################################################################
    ##########  My Mother Teaches Me Bullfighting                    ##########
    ###########################################################################

    ICP <- function (input_inci, input_cfat, input_prev, totmor, Age, bprev=0,
                     fminopt =optimset(TolX=Tol,MaxIter=MaxIt)){

      ###########################################################################
      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_mort=NULL
        input_dura=NULL
        input_rrmort=NULL
        incipoprate=NULL
        input_remi=NULL
        input_remi2=NULL
        input_prev2=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }
      ###########################################################################

      for (ag in 1:(length(Age)-1)){

        ###########################################################################
        amoeba <- function(x){

          ### intermediate variables
          intermediate_l[ag]=input_inci[ag]+x[1]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*x[1]-
            2.0*input_inci[ag]*input_cfat[ag]+x[1]**2+2.0*input_cfat[ag]*x[1]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+x[1]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*x[1]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+x[1]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          input_mort[ag+1]=(state_dead[ag+1]-state_dead[ag])/py

          (input_prev[ag+1]-0.5*(state_withcondition[ag]
                                 +state_withcondition[ag+1])/py)**2

        }


        neldermead=fminsearch(fun=amoeba, x0=neldtmp2[ag], fminopt)
        neldtmp1=neldermead$simplexopt
        neldtmp2[ag+1]=neldtmp1$x[1,]

        input_remi[1]=0

        if (neldtmp2[ag+1]>=0 & abs(neldtmp2[ag+1]-neldtmp2[ag])<0.2)   {  #
          input_remi[ag+1]=neldtmp2[ag+1]
        } else {input_remi[ag+1]=input_remi[ag]}
        ###########################################################################


        ###########################################################################
        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            intermediate_mortotaal[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {intermediate_mortotaal[ag+1]=intermediate_mortotaal[ag]}

          input_prev2[ag+1]=0.5*(state_withcondition[ag]
                                 +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }
        ###########################################################################

        intermediate_mortotaal[1]=0
        input_prev2[1]=0
      }	#For-loop end#
      remi_spline=smooth.spline(x=Age,y=input_remi,spar=0.4)
      input_remi1=predict(remi_spline,Age)
      input_remidownhill=input_remi1[[2]]
      for (ag in 1:(length(Age)-1)){
        input_remidownhill[1]=0
        if (input_remidownhill[ag+1]<0) {
          input_remidownhill[ag+1]=input_remidownhill[ag]
        }}
      output=(list(Incidence=input_inci,remissiondownhill=input_remidownhill,
                   Case_Fatality=input_cfat,Prevalence_proportion=input_prev,
                   Mortality_rate=intermediate_mortotaal,RR_Mortality=input_rrmort))
      class(output)="Epi_final"
      return(output)
    }


    ###########################################################################
    ###########################################################################
    ###########################################################################
    Final_estimates=ICP(input_inci=Compartmental_model$input_inci,
                        input_cfat=Compartmental_model$input_cfat,
                        input_prev=Compartmental_model$input_prev,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }


  if(!is.null(Compartmental_model$input_inci)&!is.null(Compartmental_model$input_cfat)&
     !is.null(Compartmental_model$input_mort)){

    ###########################################################################
    ##########  I Play Pinochle with a Horse                         ##########
    ###########################################################################

    ICM <- function (input_inci, input_cfat, input_mort, totmor, Age, bprev=0,
                     fminopt = optimset(TolX=Tol, MaxIter=MaxIt)){

      ###########################################################################
      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_prev=NULL
        input_dura=NULL
        input_rrmort=NULL
        incipoprate=NULL
        input_remi=NULL
        input_remi2=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }
      ###########################################################################

      for (ag in 1:(length(Age)-1)){

        ###########################################################################
        amoeba <- function(x){

          ### intermediate variables
          intermediate_l[ag]=input_inci[ag]+x[1]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*x[1]-
            2.0*input_inci[ag]*input_cfat[ag]+x[1]**2+2.0*input_cfat[ag]*x[1]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+x[1]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*x[1]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+x[1]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          (input_mort[ag+1]-((state_dead[ag+1]-state_dead[ag])/py))**2
        }



        neldermead=fminsearch(fun=amoeba, x0=neldtmp2[ag], fminopt)
        neldtmp1=neldermead$simplexopt
        neldtmp2[ag+1]=neldtmp1$x[1,]

        input_remi[1]=0

        if (neldtmp2[ag+1]>=0 & abs(neldtmp2[ag+1]-neldtmp2[ag])<0.2)   {  #
          input_remi[ag+1]=neldtmp2[ag+1]
        } else {input_remi[ag+1]=input_remi[ag]}
        ###########################################################################


        ###########################################################################
        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            intermediate_mortotaal[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {intermediate_mortotaal[ag+1]=intermediate_mortotaal[ag]}

          input_prev[ag+1]=0.5*(state_withcondition[ag]
                                +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }
        ###########################################################################

        intermediate_mortotaal[1]=0
        input_prev[1]=0
      }	#For-loop end#
      remi_spline=smooth.spline(x=Age,y=input_remi,spar=0.4)
      input_remi1=predict(remi_spline,Age)
      input_remidownhill=input_remi1[[2]]
      for (ag in 1:(length(Age)-1)){
        input_remidownhill[1]=0
        if (input_remidownhill[ag+1]<0) {
          input_remidownhill[ag+1]=input_remidownhill[ag]
        }}
      output=(list(Incidence=input_inci,remissiondownhill=input_remidownhill,
                   Case_Fatality=input_cfat,Prevalence_proportion=input_prev,
                   Mortality_rate=intermediate_mortotaal,RR_Mortality=input_rrmort))
      class(output)="Epi_final"
      return(output)
    }

    ###########################################################################
    ###########################################################################
    ###########################################################################
    Final_estimates=ICM(input_inci=Compartmental_model$input_inci,
                        input_cfat=Compartmental_model$input_cfat,
                        input_mort=Compartmental_model$input_mort,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }


  if(!is.null(Compartmental_model$input_remi)&!is.null(Compartmental_model$input_cfat)&
     !is.null(Compartmental_model$input_prev)){

    ###########################################################################
    ##########  I Become Supreme Lord of the Bathroom                ##########
    ###########################################################################

    RCP <- function (input_remi, input_cfat, input_prev, totmor, Age, bprev=0,
                     fminopt=optimset(TolX=Tol,MaxIter=MaxIt)){

      ###########################################################################
      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_mort=NULL
        input_dura=NULL
        input_rrmort=NULL
        incipoprate=NULL
        input_inci=NULL
        input_inci2=NULL
        input_prev2=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }
      ###########################################################################

      for (ag in 1:(length(Age)-1)){

        ###########################################################################
        amoeba <- function(x){

          ### intermediate variables
          intermediate_l[ag]=x[1]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=x[1]**2+2.0*x[1]*input_remi[ag]-
            2.0*x[1]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          input_mort[ag+1]=(state_dead[ag+1]-state_dead[ag])/py

          (input_prev[ag+1]-0.5*(state_withcondition[ag]
                                 +state_withcondition[ag+1])/py)**2

        }


        neldermead=fminsearch(fun=amoeba, x0=neldtmp2[ag], fminopt)
        neldtmp1=neldermead$simplexopt
        neldtmp2[ag+1]=neldtmp1$x[1,]

        input_inci[1]=0

        if (neldtmp2[ag+1]>=0 & abs(neldtmp2[ag+1]-neldtmp2[ag])<0.2)   {  #
          input_inci[ag+1]=neldtmp2[ag+1]
        } else {input_inci[ag+1]=input_inci[ag]}
        ###########################################################################


        ###########################################################################
        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            intermediate_mortotaal[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {intermediate_mortotaal[ag+1]=intermediate_mortotaal[ag]}

          input_prev2[ag+1]=0.5*(state_withcondition[ag]
                                 +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }
        ###########################################################################

        intermediate_mortotaal[1]=0
        input_prev2[1]=0
      }	#For-loop end#
      inci_spline=smooth.spline(x=Age,y=input_inci,spar=0.4)
      input_inci1=predict(inci_spline,Age)
      input_incidownhill=input_inci1[[2]]
      for (ag in 1:(length(Age)-1)){
        input_incidownhill[1]=0
        if (input_incidownhill[ag+1]<0) {
          input_incidownhill[ag+1]=input_incidownhill[ag]
        }}
      output=(list(incidencedownhill=input_incidownhill,Remission=input_remi,
                   Case_Fatality=input_cfat,Prevalence_proportion=input_prev,
                   Mortality_rate=intermediate_mortotaal,RR_Mortality=input_rrmort))
      class(output)="Epi_final"
      return(output)
    }


    ###########################################################################
    ###########################################################################
    ###########################################################################
    Final_estimates=RCP(input_remi=Compartmental_model$input_remi,
                        input_cfat=Compartmental_model$input_cfat,
                        input_prev=Compartmental_model$input_prev,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }


  if(!is.null(Compartmental_model$input_remi)&!is.null(Compartmental_model$input_cfat)&
     !is.null(Compartmental_model$input_mort)){

    ###########################################################################
    ##########  My Dinner Goes Up in Smoke                           ##########
    ###########################################################################

    RCM <- function (input_remi, input_cfat, input_mort, totmor, Age, bprev=0,
                     fminopt = optimset(TolX=Tol, MaxIter=MaxIt)){

      {
        ###For the Dynamic version.
        dynamic_inci=1
        dynamic_cfat=1
        dynamic_remi=1

        ###Declarations
        intermediate_l=NULL
        intermediate_q=NULL
        intermediate_w=NULL
        intermediate_v=NULL
        tmp=NULL
        state_withcondition=NULL
        state_healthy=NULL
        state_dead=NULL
        intermediate_mortotaal=NULL
        py=NULL
        input_prev=NULL
        input_dura=NULL
        input_rrmort=NULL
        incipoprate=NULL
        input_inci=NULL
        input_inci2=NULL
        neldtmp2=NULL
        neldtmp2[1]=0

        startcoh=1000.0;
        doublehigh=1.7E308;
        tiny=1.0E10-45
        state_withcondition[1]=bprev*startcoh
        state_healthy[1]=(1.0-bprev)*startcoh
        state_dead[1]=0.0

      }

      for (ag in 1:(length(Age)-1)){

        amoeba <- function(x){

          ### intermediate variables
          intermediate_l[ag]=x[1]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=x[1]**2+2.0*x[1]*input_remi[ag]-
            2.0*x[1]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          (input_mort[ag+1]-((state_dead[ag+1]-state_dead[ag])/py))**2
        }



        neldermead=fminsearch(fun=amoeba, x0=neldtmp2[ag], fminopt)
        neldtmp1=neldermead$simplexopt
        neldtmp2[ag+1]=neldtmp1$x[1,]

        input_inci[1]=0

        if (neldtmp2[ag+1]>=0 & abs(neldtmp2[ag+1]-neldtmp2[ag])<0.2)   {  #
          input_inci[ag+1]=neldtmp2[ag+1]
        } else {input_inci[ag+1]=input_inci[ag]}

        {
          ### la, qa, wa, and va
          intermediate_l[ag]=input_inci[ag]+input_remi[ag]+input_cfat[ag]
          tmp[ag]=input_inci[ag]**2+2.0*input_inci[ag]*input_remi[ag]-
            2.0*input_inci[ag]*input_cfat[ag]+input_remi[ag]**2+2.0*input_cfat[ag]*input_remi[ag]+input_cfat[ag]**2
          if (tmp[ag]>=0.0) {
            intermediate_q[ag]=sqrt(tmp[ag])
          } else {intermediate_q[ag]=0.0}
          intermediate_w[ag]=exp(-0.5*(intermediate_l[ag]+intermediate_q[ag]))
          intermediate_v[ag]=exp(-0.5*(intermediate_l[ag]-intermediate_q[ag]))



          ###nexthealthy(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_healthy[ag+1]=(2.0*(intermediate_v[ag]-intermediate_w[ag])*

                                   (state_healthy[ag]*
                                      (input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)+ #NOTE: we still need dynamic_cfat and dynamic_remi
                                      state_withcondition[ag]*input_remi[ag]*dynamic_remi)+

                                   state_healthy[ag]*
                                   (intermediate_v[ag]*(intermediate_q[ag]-intermediate_l[ag])+
                                      intermediate_w[ag]*(intermediate_q[ag]+intermediate_l[ag])))/

              (2.0*intermediate_q[ag])
          } else {state_healthy[ag+1]=state_healthy[ag]}


          ###nextdiseased(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_withcondition[ag+1]=-((intermediate_v[ag]-intermediate_w[ag])*

                                          (2.0*((input_cfat[ag]*dynamic_cfat+input_remi[ag]*dynamic_remi)*
                                                  (state_healthy[ag]+state_withcondition[ag])+
                                                  state_healthy[ag]*(-intermediate_l[ag]))
                                           -state_withcondition[ag]*intermediate_l[ag])

                                        -state_withcondition[ag]*intermediate_q[ag]*
                                          (intermediate_v[ag]+intermediate_w[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_withcondition[ag+1]=state_withcondition[ag]}



          ###nextdead(ag:integer);
          if (intermediate_q[ag]>0.0) {
            state_dead[ag+1]=((intermediate_v[ag]-intermediate_w[ag])*

                                (2.0*(state_withcondition[ag]*(input_cfat[ag]*dynamic_cfat))-
                                   intermediate_l[ag]*(state_healthy[ag]+state_withcondition[ag]))-

                                intermediate_q[ag]*(state_healthy[ag]+state_withcondition[ag])*
                                (intermediate_v[ag]+intermediate_w[ag])+

                                2.0*intermediate_q[ag]*
                                (state_healthy[ag]+state_withcondition[ag]+state_dead[ag]))/

              (2.0*intermediate_q[ag])
          } else {state_dead[ag+1]=state_dead[ag]}



          ### Prevalence and Mortality
          py=0.5*(state_healthy[ag]+state_withcondition[ag]
                  +state_healthy[ag+1]+state_withcondition[ag+1])

          if (state_dead[ag+1]-state_dead[ag]>=0) {
            intermediate_mortotaal[ag+1]=
              (state_dead[ag+1]-state_dead[ag])/py
          } else {intermediate_mortotaal[ag+1]=intermediate_mortotaal[ag]}

          input_prev[ag+1]=0.5*(state_withcondition[ag]
                                +state_withcondition[ag+1])/py

          input_rrmort[ag]=1.0+input_cfat[ag]/totmor[ag] #NOTE: TOTAL MORTALITY MUST BE GIVEN

          tmp1=input_cfat[ag]+input_remi[ag]
          if (tmp1<tiny) {
            input_dura[ag]=0.0
          } else {input_dura[ag]=1.0/tmp1}
          incipoprate[ag+1]=input_inci[ag]*0.5*
            (state_healthy[ag]+state_healthy[ag+1])/py
        }

        intermediate_mortotaal[1]=0
        input_prev[1]=0
      }	#For-loop end#
      inci_spline=smooth.spline(x=Age,y=input_inci,spar=0.4)
      input_inci1=predict(inci_spline,Age)
      input_incidownhill=input_inci1[[2]]
      for (ag in 1:(length(Age)-1)){
        input_incidownhill[1]=0
        if (input_incidownhill[ag+1]<0) {
          input_incidownhill[ag+1]=input_incidownhill[ag]
        }}
      output=(list(incidencedownhill=input_incidownhill,Remission=input_remi,
                   Case_Fatality=input_cfat,Prevalence_proportion=input_prev,
                   Mortality_rate=intermediate_mortotaal,RR_Mortality=input_rrmort,
                   Age=Age))
      class(output)="Epi_final"
      return(output)
    }

    Final_estimates=RCM(input_remi=Compartmental_model$input_remi,
                        input_cfat=Compartmental_model$input_cfat,
                        input_mort=Compartmental_model$input_mort,
                        totmor=Compartmental_model$totmor,
                        Age=Compartmental_model$Age, bprev=birth_prevalence)
  }

  plot.Epi_final <- function (fin){
    par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
    plot(fin$Prevalence_proportion,type="l",col="green",lwd=3,
         ylab="Hazard or Rate", xlab="Age", main="Output Variables",
         ylim=c(min(fin$Incidence,fin$Remission,fin$Case_Fatality,
                    fin$Prevalence_proportion,
                    fin$Mortality_rate),max(fin$Incidence,fin$Remission,fin$Case_Fatality,
                                            fin$Prevalence_proportion,
                                            fin$Mortality_rate)))
    lines(fin$Mortality_rate,col="red",lwd=3)
    lines(fin$Remission,col="yellow",lwd=3)
    lines(fin$Incidence,col="blue",lwd=3)

    if(!is.null(fin$Case_fatality_downhill)){
      par(new=T)
      plot(fin$Case_fatality_downhill,col="purple",axes=F,xlab=NA, ylab=NA, lwd=2,type="l")
      axis(side=4)}

    if(!is.null(fin$remissiondownhill)){
      par(new=T)
      plot(fin$remissiondownhill,col="yellow",axes=F,xlab=NA, ylab=NA, lwd=2,type="l")
      axis(side=4)}

    if(!is.null(fin$incidencedownhill)){
      par(new=T)
      plot(fin$incidencedownhill,col="blue",axes=F,xlab=NA, ylab=NA, lwd=2,type="l")
      axis(side=4)}

    legend("topright", inset=c(-0.45,0), legend=c("Incidence", "Remission", "Case Fatality", "Mortality", "Prevalence"),
           col=c("blue","yellow","purple","red", "green"), lwd=c(3:3),title="Legend:", bty="n")
  }

  plot(Final_estimates)
  plot_Final_estimates=recordPlot()

  Tabulated <- function (AgeLB, AgeUB, Final){

    Incidence1=NULL
    Remission1=NULL
    Case_Fatality1=NULL
    Mortality1=NULL
    Prevalence1=NULL
    RR_Mortality1=NULL

    for (j in 1:length(AgeUB)){
      inci=0; remi=0; cfat=0; mort=0; prev=0; rr_mort=0; count=0;
      for (i in 1:length(Final$Incidence)){
        if(i>=AgeLB[j]&i<=AgeUB[j]){
          inci=inci+Final$Incidence[i]
          remi=remi+Final$Remission[i]
          cfat=cfat+Final$Case_fatality_downhill[i]
          mort=mort+Final$Mortality_rate[i]
          prev=prev+Final$Prevalence_proportion[i]
          rr_mort=rr_mort+Final$RR_Mortality[i]
          count=count+1
        }}
      Incidence1[j]=round(inci/count,digits=4)
      Remission1[j]=round(remi/count,digits=4)
      Case_Fatality1[j]=round(cfat/count,digits=4)
      Mortality1[j]=round(mort/count,digits=4)
      Prevalence1[j]=round(prev/count,digits=4)
      RR_Mortality1[j]=round(rr_mort/count,digits=4)
    }
    Tabulated = list(AgeLB=AgeLB, AgeUB=AgeUB, Incidence=Incidence1, Remission=Remission1,
                     Case_Fatality=Case_Fatality1,Mortality=Mortality1, Prevalence=Prevalence1,
                     RR_Mortality=RR_Mortality1)
    return(Tabulated)
  }

  Final_Table=Tabulated(AgeLB=Data$AgeLB, AgeUB=Data$AgeUB, Final=Final_estimates)

  Epi_data3 <- function (Output){

    Age1=NULL
    Age2=NULL
    Incidence1=NULL
    Incidence2=NULL
    Remission1=NULL
    Remission2=NULL
    Case_Fatality1=NULL
    Case_Fatality2=NULL
    Mortality1=NULL
    Mortality2=NULL
    Prevalence1=NULL
    Prevalence2=NULL
    RR_Mortality1=NULL
    RR_Mortality2=NULL
    for (i in 1:length(Output$AgeUB)){
      for (j in 1:Output$AgeUB[i]){
        if(j>=Output$AgeLB[i]){
          Age1[j]=j
          Incidence1[j]=Output$Incidence[i]
          Remission1[j]=Output$Remission[i]
          Case_Fatality1[j]=Output$Case_Fatality[i]
          Mortality1[j]=Output$Mortality[i]
          Prevalence1[j]=Output$Prevalence[i]
          RR_Mortality1[j]=Output$RR_Mortality[i]
        }}}
    for (ag in 1:length(Age1)){
      Age2[ag+1]=Age1[ag]
      Incidence2[ag+1]=Incidence1[ag]
      Remission2[ag+1]=Remission1[ag]
      Case_Fatality2[ag+1]=Case_Fatality1[ag]
      Mortality2[ag+1]=Mortality1[ag]
      Prevalence2[ag+1]=Prevalence1[ag]
      RR_Mortality2[ag+1]=RR_Mortality1[ag]
    }

    Age2[1]=0
    if (!is.null(Incidence2)){Incidence2[1]=0}
    if (!is.null(Remission2)){Remission2[1]=0}
    if (!is.null(Case_Fatality2)){Case_Fatality2[1]=0}
    if (!is.null(Mortality2)){Mortality2[1]=0}
    if (!is.null(Prevalence2)){Prevalence2[1]=0}
    if (!is.null(RR_Mortality2)){RR_Mortality2[1]=0}
    output=list(Age=Age2,Incidence=Incidence2,Remission=Remission2,
                Case_Fatality=Case_Fatality2,Mortality=Mortality2,
                Prevalence=Prevalence2, RR_Mortality=RR_Mortality2,
                Total_Mortality=Total_Mortality)
    class(output) <- "Epi_data3"
    return(output)
  }

  plot.Epi_data3 <- function (Epi_raw){
    par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
    plot(Epi_raw$Age,Epi_raw$Age,type="n",col=4,lwd=3,
         ylim=c(min(c(Epi_raw$Incidence,Epi_raw$Remission,Epi_raw$Mortality,
                      Epi_raw$Case_Fatality,Epi_raw$Prevalence)),
                max(c(Epi_raw$Incidence,Epi_raw$Remission,Epi_raw$Mortality,
                      Epi_raw$Case_Fatality,Epi_raw$Prevalence))),
         ylab=c("Hazard or Rate"),xlab="Age",main="Epidemiology of a Disease")
    count=0
    legend("topright", inset=c(-0.2,count), legend=c(""), title="Legend:", bty="n")
    if (length(Epi_raw$Age)==length(Epi_raw$Incidence)){
      lines(Epi_raw$Age,Epi_raw$Incidence,col="blue",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.3,count),legend=c("Incidence"),lwd=3, col="blue", bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Remission)){
      lines(Epi_raw$Age,Epi_raw$Remission,col="yellow",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.32,count),legend=c("Remission"),lwd=3,col="yellow",bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Mortality)){
      lines(Epi_raw$Age,Epi_raw$Mortality,col="red",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.285,count),legend=c("Mortality"),lwd=3,col="red",bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Case_Fatality)){
      lines(Epi_raw$Age,Epi_raw$Case_Fatality,col="purple",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.358,count),legend=c("Case Fatality"),lwd=3,col="purple",bty="n")}
    if (length(Epi_raw$Age)==length(Epi_raw$Prevalence)){
      lines(Epi_raw$Age,Epi_raw$Prevalence,col="green",lwd=3)
      count=count+0.1
      legend("topright",inset=c(-0.33,count),legend=c("Prevalence"),lwd=3,col="green",bty="n")}
  }

  Epi_raw3=Epi_data3(Output=Final_Table)
  plot(Epi_raw3)
  plot_Final_estimates_group=recordPlot()

  Output=list("Input"=Input,"Plot_input"=plot_input, "Cubic_spline"=Cubic_spline_data,
              "Plot_cubic_spline"=plot_cubic_spline, "Final_estimates"=Final_estimates,
              "Plot_Final_estimates"=plot_Final_estimates, "Final_Table"=Final_Table,
              "Plot_Final_groupestimates"=plot_Final_estimates_group)
  return(Output)
}


