!######################################################################
!# MODULE: mcmcdata
!#
!# AUTHORS:
!#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
!#         Rob Deardon <robert.deardon@ucalgary.ca>
!#
!# DESCRIPTION:
!#
!#     Markovchain Monte Carlo (MCMC) algorithm for the estimation of
!#     individual-level model (ILM) parameters with two disease types: SI and SIR
!#
!# Algorithm based on:
!#
!#     Deardon R, Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley,
!#     M. J., Savill, N. J., Shaw, D. J.,  Woolhouse, M. E. (2010).
!#     Inference for individual level models of infectious diseases in large
!#     populations. Statistica Sinica, 20, 239-261.
!#
!#     This program is free software; you can redistribute it and/or
!#     modify it under the terms of the GNU General Public License,
!#     version 2, as published by the Free Software Foundation.
!#
!#     This program is distributed in the hope that it will be useful,
!#     but without any warranty; without even the implied warranty of
!#     merchantability or fitness for a particular purpose. See the GNU
!#     General Public License, version 2, for more details.
!#
!#     A copy of the GNU General Public License, version 2, is available
!#     at http://www.r-project.org/Licenses/GPL-2
!#
!# Part of the R/EpiILM package
!# Contains:
!#            mcmc           .................. subroutine
!#            conmcmc        .................. subroutine
!#            initrandomseed .................. subroutine
!#            half_normal    .................. function
!#            gamma_density  .................. function
!#            rand_normal    .................. function
!#            like           .................. subroutine
!#            likesir        .................. subroutine
!#            likecon        .................. subroutine
!#            likeconsir     .................. subroutine
!######################################################################

    module mcmcdata
    use ISO_C_BINDING
    implicit none
    public:: mcmc, conmcmc

    contains

!######################################################################

    subroutine mcmc(tnum, x, y, tau, n, lambda, tmin, tmax, nsim, aalpha, ns, ni, &
        & bbeta, covmat, prostda, prostdb, anum, bnum, halfvar, unifmin, unifmax, &
        & gshape, gscale, halfvarb, unifminb, unifmaxb, gshapeb, gscaleb, simalpha, &
        & simbeta, sspark, flag, prostdsp, snum, halfvarsp, unifminsp, unifmaxsp, &
        & gshapesp, gscalesp, simspark, llikeval, tempseed) bind(C, name="mcmc_")
    !MCMC algorithm for parameter estimation of purely spatial ILMs: SI and SIR
    implicit none

    !Declarations
    integer (C_INT), intent(in)  :: n, nsim, tmax, ns, ni, tmin, tempseed
    integer (C_INT), intent(in)  :: tnum, anum, bnum      !likelihood, model and prior selections
    integer (C_INT), intent(in)  :: tau(n), lambda(n)     !infection times and infperiod
    real (C_DOUBLE), intent(in)  :: aalpha(ns), bbeta(ni) !initial values of parameters
    real (C_DOUBLE), intent(in)  :: prostda(ns), prostdb  !proposal std. deviation
    real (C_DOUBLE), intent(in)  :: halfvar(ns), gshape(ns), gscale(ns) !prior parameters
    real (C_DOUBLE), intent(in)  :: unifmin(ns), unifmax(ns)            !prior parameters
    real (C_DOUBLE), intent(in)  :: halfvarb, gshapeb, gscaleb, unifminb, unifmaxb !prior parameters
    real (C_DOUBLE), intent(in)  :: x(n), y(n), covmat(n,ns)  !locations and covariates
    real (C_DOUBLE), intent(out) :: simalpha(nsim,ns), simbeta(nsim,ni), simspark(nsim) !result
    real (C_DOUBLE), intent(out) :: llikeval(nsim)                                      !result
    integer (C_INT), intent(in)  :: flag, snum                     !flag for spark term
    real (C_DOUBLE), intent(in)  :: sspark, prostdsp               !initial value and proposal std. dev
    real (C_DOUBLE), intent(in)  :: halfvarsp, unifminsp           !prior parameters for spark
    real (C_DOUBLE), intent(in)  :: unifmaxsp, gshapesp, gscalesp  !prior parameters for spark

    double precision :: alpha(ns), beta(ni)
    double precision :: z_alpha, z_beta, alpha_n(ns), beta_n(ni)
    double precision :: ratio, ratio_b, value, value_ini
    double precision :: u, psi, psib
    double precision :: value_b,  value_fb, value_fb_ini
    double precision :: value_f, value_f_ini
    integer          :: m, i
    double precision :: spark, z_sp, spark_n, ratio_sp, psisp
    double precision :: value_sp, value_fsp, value_fsp_ini

    !initialzing random seed
    if (tempseed .NE. 0) then
        call initrandomseed(tempseed)
    end if

    !initializations for acceptance probability calculations
    psi   = 0.d0
    psib  = 0.0d0
    psisp = 0.0d0

    !initial value assignments
    do i = 1, ns
      alpha(i) = aalpha(i)
    end do
    beta(ni) = bbeta(ni)
    spark    = sspark

    do i = 1,ns
      simalpha(1,i) = alpha(i)
    end do
    simbeta(1,ni) = beta(ni)
    simspark(1)   = spark

    do i = 1, ns
      alpha_n(i) = alpha(i)
    end do

    !log-likelihood value selection with initial parameter values
    SELECT CASE (tnum)
    CASE (1) !SI model
    call like(x, y, tau, n, tmin, tmax, ns, ni, alpha, beta, spark, covmat, value_ini)

    CASE (2) !SIR model
    call likesir(x, y, tau, lambda, n, tmin, tmax, ns, ni, alpha, beta, spark, covmat, value_ini)
    END SELECT

    llikeval(1) = value_ini   !initial log-likelihood value

    !simulation starts here
    do m = 1, (nsim-1)
      !estimation of alpha parameter (s)
      do i = 1, ns                                  !do-loop for each alpha

        z_alpha    = rand_normal(0.0d0, prostda(i)) !generate random numbers from proposal
        alpha_n(i) = alpha(i) + z_alpha             !update alpha parameter (s)

        if (alpha_n(i) .GT. 0.0d0) then             !model condition for alpha parameter (s)

          SELECT CASE (tnum)
          CASE (1) !SI log-likelihood
          call like(x, y, tau, n, tmin, tmax, ns, ni, alpha_n, beta, spark, covmat, value)

          CASE (2) !SIR log-likelihood
          call likesir(x, y, tau, lambda, n, tmin, tmax, ns, ni, alpha_n, beta, spark, covmat, value)
          END SELECT

          !prior selections for alpha parameter estimation
          SELECT CASE (anum)
          CASE (1) !gamma prior
          value_f     = gamma_density(alpha_n(i), gshape(i), gscale(i))
          value_f_ini = gamma_density(alpha(i), gshape(i), gscale(i))

          !calcualte acceptance probability
          ratio = (value - value_ini) + (value_f - value_f_ini)
          psi   = min(1.0d0, exp(ratio))

          CASE (2) !half normal prior
          value_f     = half_normal(alpha_n(i), halfvar(i))
          value_f_ini = half_normal(alpha(i), halfvar(i))

          !calcualte acceptance probability
          ratio = (value - value_ini) + (value_f - value_f_ini)
          psi   = min(1.0d0, exp(ratio))

          CASE (3) !uniform prior
          if ((unifmin(i) .LT. alpha_n(i)) .AND. (alpha_n(i).LT. unifmax(i))) then

            !calcualte acceptance probability
            ratio = (value - value_ini)
            psi   = min(1.0d0, exp(ratio))
          else
            psi = -1.0d0
          end if
          END SELECT

          !acceptance probbaility check for alpha parameters
          call random_number(u)
          if (psi >= u) then
            simalpha(m+1,i) = alpha_n(i)
            value_ini = value
          else
            simalpha(m+1,i) = alpha(i)
          end if
        else
          simalpha(m+1,i) = alpha(i)
        end if                        !model condition for alpha parameter (s) ends here
        alpha_n(i) = simalpha(m+1,i)  !updating initial value
     end do                          !do-loop for each alpha ends here

      !estimation of beta parameter
      z_beta     = rand_normal(0.0d0, prostdb) !generate from proposal
      beta_n(ni) = beta(ni) + z_beta           !update beta

      if (beta_n(ni) > 0.0d0) then             !model condition for beta parameter

        SELECT CASE (tnum)
        CASE (1) !SI log-likelihood
        call like(x, y, tau, n, tmin, tmax, ns, ni, alpha_n, beta_n, spark, covmat, value_b)

        CASE (2) !SIR log-likelihood
        call likesir(x, y, tau, lambda, n, tmin, tmax, ns, ni, alpha_n, beta_n, spark, covmat, value_b)

        END SELECT

        !prior selections for beta parameter estimation
        SELECT CASE (bnum)
        CASE (1)  !gamma prior
        value_fb     = gamma_density(beta_n(ni), gshapeb, gscaleb)
        value_fb_ini = gamma_density(beta(ni), gshapeb, gscaleb)

        !calcualte acceptance probability
        ratio_b = (value_b - value_ini) + (value_fb - value_fb_ini)
        psib    = min(1.0d0, exp(ratio_b))

        CASE (2)  !half normal prior
        value_fb     = half_normal(beta_n(ni), halfvarb)
        value_fb_ini = half_normal(beta(ni), halfvarb)

        !calcualte acceptance probability
        ratio_b = (value_b - value_ini) + (value_fb - value_fb_ini)
        psib    = min(1.0d0, exp(ratio_b))

        CASE (3) !uniform prior
        if ((unifminb .LT. beta_n(ni)).AND. (beta_n(ni) .LT. unifmaxb)) then

          !calcualte acceptance probability
          ratio_b = value_b - value_ini
          psib    = min(1.0d0, exp(ratio_b))
        else
          psib = -1.0d0
        end if
        END SELECT

        !acceptance probbaility check for beta parameter
        call random_number(u)
        if (psib >= u) then
          simbeta(m+1,ni) = beta_n(ni)
          value_ini = value_b
        else
          simbeta(m+1,ni) = beta(ni)
        end if
      else
        simbeta(m+1,ni) = beta(ni)
      end if                       !model condition for beta parameter ends here
      beta_n(ni) = simbeta(m+1,ni) !update initial value

      !estimation of spark parameter
      if (flag .EQ. 1) then                    !if spark parameter exists in the model
        z_sp    = rand_normal(0.0d0, prostdsp) !generate from proposal
        spark_n = spark + z_sp                 !update spark parameter

        if (spark_n > 0.0d0) then              !model condition for spark parameter

          SELECT CASE (tnum)
          CASE (1) !SI log-likelihood
          call like(x, y, tau, n, tmin, tmax, ns, ni, alpha_n, beta_n, spark_n, covmat, value_sp)

          CASE (2) !SIR log-likelihood
          call likesir(x, y, tau, lambda, n, tmin, tmax, ns, ni, alpha_n, beta_n, spark_n, covmat, value_sp)

          END SELECT

          !prior distribution selection for spark parameter
          SELECT CASE (snum)
          CASE (1)  !gamma prior
          value_fsp     = gamma_density(spark_n, gshapesp, gscalesp)
          value_fsp_ini = gamma_density(spark, gshapesp, gscalesp)

          !calcualte acceptance probability
          ratio_sp = (value_sp - value_ini) + (value_fsp - value_fsp_ini)
          psisp    = min(1.0d0, exp(ratio_sp))

          CASE (2)  !half normal prior
          value_fsp     = half_normal(spark_n, halfvarsp)
          value_fsp_ini = half_normal(spark, halfvarsp)

          !calcualte acceptance probability
          ratio_sp = (value_sp - value_ini) + (value_fsp - value_fsp_ini)
          psisp    = min(1.0d0, exp(ratio_sp))

          CASE (3)  !uniform prior
          if ((unifminsp .LT. spark_n).AND. (spark_n .LT. unifmaxsp)) then

            !calcualte acceptance probability
            ratio_sp = value_sp - value_ini
            psisp    = min(1.0d0, exp(ratio_sp))
          else
            psisp = -1.0d0
          end if
          END SELECT

          !acceptance probbaility check for spark
          call random_number(u)
          if (psisp >= u) then
            simspark(m+1) = spark_n
          else
            simspark(m+1) = spark
          end if
        else
          simspark(m+1) = spark
        end if                     !model condition for spark parameter ends here
          spark_n = simspark(m+1)  !update initial
      else
        spark   = 0.0d0
        spark_n = 0.0d0
      end if                       !flag ends here

      SELECT CASE (tnum)
      CASE (1) !SI log-likelihood
      call like(x, y, tau, n, tmin, tmax, ns, ni, alpha_n, beta_n, spark_n, covmat, value_ini)
      llikeval(m+1) = value_ini

      CASE (2)  !SIR log-likelihood
      call likesir(x, y, tau, lambda, n, tmin, tmax, ns, ni, alpha_n, beta_n, spark_n, covmat, value_ini)
      llikeval(m+1) = value_ini
      END SELECT

      do i = 1, ns
       alpha(i) = alpha_n(i)
      end do
      beta(ni) = beta_n(ni)
      spark    = spark_n
    end do               !simulation ends
    end subroutine mcmc  !subroutine ends

!######################################################################

    subroutine conmcmc(tnum, tau, n, lambda, tmin, tmax, nsim, aalpha, ns, ni, &
   		              & bbeta, covmat, network, prostda, prostdb, anum, bnum, halfvar, unifmin, &
                      & unifmax, gshape, gscale, halfvarb, unifminb, unifmaxb, gshapeb, gscaleb, &
                      & simalpha, simbeta, sspark, flag, prostdsp, snum, halfvarsp, &
                      & unifminsp, unifmaxsp, gshapesp, gscalesp, simspark, &
                      & llikeval, tempseed) bind(C, name="conmcmc_")
    !MCMC algorithm for parameter estimation of contact network based ILMs: SI and SIR
    implicit none

    !Declarations
    integer (C_INT), intent(in)  :: n, nsim, tmax, ns, ni, tmin, tempseed
    integer (C_INT), intent(in)  :: tnum, anum, bnum          !likelihood, model and prior selections
    integer (C_INT), intent(in)  :: lambda(n), tau(n)         !infperiod and infection times
    real (C_DOUBLE), intent(in)  :: aalpha(ns), bbeta(ni)     !initial parameter values
    real (C_DOUBLE), intent(in)  :: prostda(ns), prostdb(ni)  !proposal std. deviation
    real (C_DOUBLE), intent(in)  :: halfvar(ns), gshape(ns), gscale(ns)    !prior parameter values
    real (C_DOUBLE), intent(in)  :: unifmin(ns), unifmax(ns)               !prior parameter values
    real (C_DOUBLE), intent(in)  :: halfvarb(ni), gshapeb(ni), gscaleb(ni) !prior parameter values
    real (C_DOUBLE), intent(in)  :: unifminb(ni), unifmaxb(ni)             !prior parameter values
    real (C_DOUBLE), intent(in)  :: network(n, n, ni)                      !contact network
    real (C_DOUBLE), intent(in)  :: covmat(n, ns)                          !covariates
    real (C_DOUBLE), intent(out) :: simalpha(nsim,ns), simbeta(nsim,ni), simspark(nsim) !result
    real (C_DOUBLE), intent(out) :: llikeval(nsim)                                      !result

    integer (C_INT), intent(in) :: flag, snum                      !flag for spark term
    real (C_DOUBLE), intent(in) :: sspark, prostdsp                !initial and proposal std. dev
    real (C_DOUBLE), intent(in) :: halfvarsp, unifminsp, unifmaxsp !prior parameters for spark
    real (C_DOUBLE), intent(in) :: gshapesp, gscalesp              !prior parameters for spark

    double precision :: alpha(ns), beta(ni)
    double precision :: z_alpha, z_beta, alpha_n(ns), beta_n(ni)
    double precision :: ratio, ratio_b, value, value_ini
    double precision :: u, psi, psib
    double precision :: value_b, value_fb, value_fb_ini
    double precision :: value_f, value_f_ini
    integer          :: m, i
    double precision :: spark,z_sp, spark_n, ratio_sp, psisp
    double precision :: value_sp, value_fsp, value_fsp_ini

    !initialzing random seed
    if (tempseed .NE. 0) then
        call initrandomseed(tempseed)
    end if

    !initializations for acceptance probability calculations
    psi   = 0.d0
    psib  = 0.0d0
    psisp = 0.0d0

    !initial value assignment
    do i = 1, ns
      alpha(i)      = aalpha(i)
      simalpha(1,i) = alpha(i)
      alpha_n(i)    = alpha(i)
    end do

    do i = 1, ni
      beta(i)      = bbeta(i)
      simbeta(1,i) = beta(i)
      beta_n(i)    = beta(i)
    end do

    spark       = sspark
    simspark(1) = spark

    !log-likelihood value selection with initial parameter values
    SELECT CASE (tnum)
    CASE (1) !SI model
    call likecon(tau, n, ns, ni, tmin, tmax, alpha, beta, spark, covmat, network, value_ini)

    CASE (2)  !SIR model
    call likeconsir(tau, lambda, n, ns, ni, tmin, tmax, alpha, beta, spark, covmat, network, value_ini)
    END SELECT

    llikeval(1) = value_ini  !initial log-likelihood

    !simulation starts here
    do m = 1, (nsim-1)
      !estimation of alpha parameter (s)
      do i = 1, ns                                  !do-loop for each alpha
        z_alpha    = rand_normal(0.0d0, prostda(i)) !generate random numbers from proposal
        alpha_n(i) = alpha(i) + z_alpha             !update alpha

        if (alpha_n(i) .GT. 0.0d0) then             !model condition for alpha parameter(s)

        !log-likelihood value selection
          SELECT CASE (tnum)
          CASE (1) !SI model
          call likecon(tau, n, ns, ni, tmin, tmax, alpha_n, beta, spark, covmat, network, value)

          CASE (2)  !SIR model
          call likeconsir(tau, lambda, n, ns, ni, tmin, tmax, alpha_n, beta, spark, covmat, network, value)
          END SELECT

          !prior distribution selection
          SELECT CASE (anum)
          CASE (1)  !gamma prior
          value_f     = gamma_density(alpha_n(i), gshape(i), gscale(i))
          value_f_ini = gamma_density(alpha(i), gshape(i), gscale(i))

          !calcualte acceptance probbaility
          ratio = (value - value_ini) + (value_f - value_f_ini)
          psi   = min(1.0d0, exp(ratio))

          CASE (2) !half normal prior
          value_f     = half_normal(alpha_n(i), halfvar(i))
          value_f_ini = half_normal(alpha(i), halfvar(i))

          !calcualte acceptance probbaility
          ratio = (value - value_ini) + (value_f - value_f_ini)
          psi   = min(1.0d0, exp(ratio))

          CASE (3)  !uniform prior
          if ((unifmin(i) .LT. alpha_n(i)).AND. (alpha_n(i) .LT. unifmax(i))) then

            !calcualte acceptance probbaility
            ratio = (value - value_ini)
            psi   = min(1.0d0, exp(ratio))
          else
            psi = -1.0d0
          end if
          END SELECT

          !acceptance probbaility check for alpha
          call random_number(u)
          if (psi >= u) then
            simalpha(m+1,i) = alpha_n(i)
            value_ini = value
          else
            simalpha(m+1,i) = alpha(i)
          end if
        else
        simalpha(m+1,i) = alpha(i)
        end if                       !model condition for alpha parameter (s) ends here

        alpha_n(i) = simalpha(m+1,i) !updating inital value
      end do                         !do-loop for each alpha ends here

      !estimation of beta parameter (s)
      do i = 1, ni
        z_beta    = rand_normal(0.0d0, prostdb(i)) !generate from proposal
        beta_n(i) = beta(i) + z_beta               !update beta

        if (beta_n(i) > 0.0d0) then                !model condition for beta parameter

          SELECT CASE (tnum)
          CASE (1)  !SI log-likelihood
          call likecon(tau, n, ns, ni, tmin, tmax, alpha_n, beta_n, spark, covmat, network, value_b)

          CASE (2)  !SIR log-likelihood
          call likeconsir(tau, lambda, n, ns, ni, tmin, tmax, alpha_n, beta_n, spark, covmat, network, value_b)
          END SELECT

          !prior distribution selection
          SELECT CASE (bnum)
          CASE (1) !gamma prior
          value_fb     = gamma_density(beta_n(i), gshapeb(i), gscaleb(i))
          value_fb_ini = gamma_density(beta(i), gshapeb(i), gscaleb(i))

          !calcualte acceptance probbaility
          ratio_b = (value_b - value_ini) + (value_fb - value_fb_ini)
          psib    = min(1.0d0, exp(ratio_b))

          CASE (2) !half normal prior
          value_fb     = half_normal(beta_n(i), halfvarb(i))
          value_fb_ini = half_normal(beta(i), halfvarb(i))

          !calcualte acceptance probbaility
          ratio_b = (value_b - value_ini) + (value_fb - value_fb_ini)
          psib    = min(1.0d0, exp(ratio_b))

          CASE (3)  !uniform prior
          if ((unifminb(i) .LT. beta_n(i)).AND. (beta_n(i).LT. unifmaxb(i))) then

            !calcualte acceptance probbaility
            ratio_b = value_b - value_ini
            psib    = min(1.0d0, exp(ratio_b))
          else
            psib = -1.0d0
          end if
          END SELECT

          !acceptance probbaility check for beta
          call random_number(u)
          if (psib >= u) then
            simbeta(m+1,i) = beta_n(i)
            value_ini = value_b
          else
            simbeta(m+1,i) = beta(i)
          end if
        else
          simbeta(m+1,i) = beta(i)
        end if                      !model condition for beta parameter (s) ends here

        beta_n(i) = simbeta(m+1,i)  !updating initial value(s)
      end do                        !do-loop for each beta ends here

      !estimation of spark parameter
      if (flag .EQ. 1) then                       !if spark parameter exists in the model

        z_sp    = rand_normal(0.0d0, prostdsp)    !generate from proposal
        spark_n = spark + z_sp                    !update spark

        if (spark_n > 0.0d0) then                 !model condition for spark

          SELECT CASE (tnum)
          CASE (1)  !SI log-likelihood
          call likecon(tau, n, ns, ni, tmin, tmax, alpha_n, beta_n, spark_n, covmat, network, value_sp)

          CASE (2)  !SIR log-likelihood
          call likeconsir(tau, lambda, n, ns, ni, tmin, tmax, alpha_n, beta_n, spark_n, covmat, network, value_sp)
          END SELECT

          !prior distribution selection
          SELECT CASE (snum)
          CASE (1) !gamma prior
          value_fsp     = gamma_density(spark_n, gshapesp, gscalesp)
          value_fsp_ini = gamma_density(spark, gshapesp, gscalesp)

          !calcualte acceptance probbaility
          ratio_sp = (value_sp - value_ini) + (value_fsp - value_fsp_ini)
          psisp    = min(1.0d0, exp(ratio_sp))

          CASE (2) !half normal prior
          value_fsp     = half_normal(spark_n, halfvarsp)
          value_fsp_ini = half_normal(spark, halfvarsp)

          !calcualte acceptance probbaility
          ratio_sp = (value_sp - value_ini) + (value_fsp - value_fsp_ini)
          psisp    = min(1.0d0, exp(ratio_sp))

          CASE (3)  !uniform prior
          if ((unifminsp .LT. spark_n).AND. (spark_n .LT. unifmaxsp)) then

            !calcualte acceptance probbaility
            ratio_sp = value_sp - value_ini
            psisp    = min(1.0d0, exp(ratio_sp))
          else
            psisp = -1.0d0
          end if
          END SELECT

         !acceptance probbaility check for spark
          call random_number(u)
          if (psisp >= u) then
            simspark(m+1) = spark_n
          else
            simspark(m+1) = spark
          end if
        else
          simspark(m+1) = spark
        end if                   !model condition for spark parameter ends here

        spark_n = simspark(m+1)  !updating initial  value
      else
        spark   = 0.0d0
        spark_n = 0.0d0
      end if                     !flag ends here

      SELECT CASE (tnum)
      CASE (1)  !SI log-likelihood
      call likecon(tau, n, ns, ni, tmin, tmax, alpha_n, beta_n, spark_n, covmat, network, value_ini)
      llikeval(m+1) = value_ini

      CASE (2)  !SIR log-likelihood
      call likeconsir(tau, lambda, n, ns, ni, tmin, tmax, alpha_n, beta_n, spark_n, covmat, network, value_ini)
      llikeval(m+1) = value_ini
      END SELECT

      do i = 1, ns
        alpha(i) = alpha_n(i)
      end do
      do i = 1, ni
        beta(i) = beta_n(i)
      end do
      spark = spark_n

     end do                   !simulation ends here
     end subroutine conmcmc   !subroutine ends

!######################################################################

    subroutine initrandomseed(tempseed)
    !The seed for the random number generation method random_number() has been reset
    implicit none

    integer :: n
    integer, intent(in) :: tempseed
    integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))
        seed = tempseed
        call random_seed(put = seed)

    end subroutine initrandomseed

!######################################################################

    FUNCTION half_normal(alpha, b) RESULT(pdf)
    !half normal density calculation
    implicit none

    double precision, parameter:: pi = 3.141592653589793D+00
    double precision :: alpha, b
    double precision :: val, pdf

    val =  sqrt (2 / pi * b) * exp ( - 0.5D+00 * ((alpha)**2 / b) )
    pdf = log(val)

    END FUNCTION  ! returns log value

!######################################################################

    FUNCTION gamma_density(x, a, b) RESULT(pdf)
    !gamma density calculation
    implicit none

    double precision :: x, a, b
    double precision :: dn, pdf

    dn = (x**(a-1)) * exp(- (x * b))
    pdf = log(dn)

    END FUNCTION ! returns log value

!######################################################################

    FUNCTION rand_normal(mean, stdev) RESULT(c)
    !generate random variable from normal distribution
    implicit none

    double precision            :: mean, stdev, c, temp(2), r, theta
    double precision, parameter :: pi = 3.141592653589793D+00

    call random_number(temp)
    r = (-2.0d0 * log(temp(1)))**0.5
    theta = 2.0d0 * pi * temp(2)
    c = mean + stdev * r * sin(theta)

    END FUNCTION

!######################################################################

    subroutine like(x, y, tau, n, tmin, tmax, ns, ni, alpha, beta, &
                    & spark, covmat, val)
    !log-likelihood for purely spatial model: SI
    implicit none

    !Declarations
    integer, intent(in)          :: n, tau(n), tmax, ns, ni, tmin
    double precision, intent(in) :: alpha(ns), beta(ni), spark  !parameters
    double precision, intent(in) :: x(n), y(n)                  !locations
    double precision, intent(in) :: covmat(n, ns)               !covariates
    double precision, intent(out):: val                         !result

    integer          :: i, j, t
    double precision :: eu(n, n), Somega(n)
    double precision :: dx1, dx2, p1, p2

    !Calculate the distance matrix
    do i = 1, n
      do j = i, n
        eu(i,j) = sqrt(((x(i)-x(j))**2) + ((y(i)-y(j))**2))
        eu(j,i) = eu(i,j)
      end do
    end do

    Somega = matmul(covmat,alpha) !susceptibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
        !infectious period
        if (tau(i)==(t+1)) then
          dx1 = 0.0d0
          do j = 1, n
            if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0)) then
              dx1 = dx1 + (eu(i,j)**(-beta(ni)))
            end if
         end do
         p1 = 1.0d0 - exp(-(Somega(i) * dx1 + spark))
         val = val + log(p1)
       end if
      !susceptible period
      if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
        dx2 = 0.0d0
        do j = 1, n
          if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0))then
            dx2 = dx2 + (eu(i,j)**(-beta(ni)))
          end if
        end do
        p2 = exp(-(Somega(i) * dx2 + spark))
        val = val + log(p2)                   !result
        end if
      end do
    end do

    end subroutine like

!######################################################################

    subroutine likesir(x, y, tau, lambda, n, tmin, tmax, ns, ni, alpha, &
                       & beta, spark, covmat, val)
    !log-likelihood for purely spatial model: SIR
    implicit none

    !Declarations
    integer, intent(in)          :: n, tau(n), lambda(n), tmax, ns, ni, tmin
    double precision, intent(in) :: alpha(ns), beta(ni), spark  !parameters
    double precision, intent(in) :: x(n), y(n)                  !locations
    double precision, intent(in) :: covmat(n, ns)               !covariates
    double precision, intent(out):: val                         !result

    integer          :: i, j, t
    double precision :: eu(n,n), Somega(n)
    double precision :: dx1, dx2, p1, p2

    !Calculate the distance matrix
    do i = 1, n
      do j = i, n
        eu(i,j) = sqrt(((x(i)-x(j))**2) + ((y(i)-y(j))**2))
        eu(j,i) = eu(i,j)
      end do
    end do

    Somega = matmul(covmat, alpha) !susceptibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
        !infectious period
        if (tau(i)==(t+1)) then
          dx1 = 0.0d0
          do j = 1, n
            if (tau(j) .NE. 0) then
              if ((tau(j) .LT. (t+1)) .AND. (tau(j) + lambda(j) .GE. (t+1))) then
                dx1 = dx1 + (eu(i,j)**(-beta(ni)))
              end if
            end if
          end do
          p1 = 1.0d0 - exp(-(Somega(i) * dx1 + spark))
          val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
          dx2 = 0.0d0
          do j = 1, n
            if (tau(j) .NE. 0) then
              if ((tau(j) .LT. (t+1)) .AND. (tau(j)+lambda(j) .GE. (t+1)))then
                dx2 = dx2 + (eu(i,j)**(-beta(ni)))
              end if
            end if
          end do
          p2 = exp(-(Somega(i) * dx2 + spark))
          val = val + log(p2)  !result
        end if
      end do
    end do

    end subroutine likesir

!######################################################################

    subroutine likecon(tau, n, ns, ni, tmin, tmax, alpha, beta, &
                      & spark, covmat, network, val)
    !log-likelihood for contact network model: SI
    implicit none

    !Declarations
    integer, intent(in)           :: n, tmax, ns, ni, tmin
    integer, intent(in)           :: tau(n)                     !infection times
    double precision, intent(in)  :: alpha(ns), beta(ni), spark !parameters
    double precision, intent(in)  :: covmat(n,ns)               !covariates
    double precision, intent(in)  :: network(n,n,ni)            !contact network
    double precision, intent(out) :: val                        !result

    integer          :: i, j, t, k
    double precision :: dx1, dx2, p1, p2, Somega(n)

    Somega = matmul(covmat,alpha) !susceptibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
        !infectious period
        if (tau(i)==(t+1)) then
          dx1 = 0.0d0
          do j = 1, n
            if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0)) then
              do k = 1, ni
                dx1 = dx1 + beta(k) * network(i,j,k)
              end do
            end if
          end do
          p1 = 1.0d0 - exp(-(Somega(i) * dx1 + spark))
          val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
          dx2 = 0.0d0
          do j = 1, n
            if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0))then
              do k = 1, ni
                dx2 = dx2 + beta(k) * network(i,j,k)
              end do
            end if
          end do
          p2 = exp(-(Somega(i) * dx2 + spark))
          val = val + log(p2)  !result
        end if
      end do
    end do

    end subroutine likecon

!######################################################################

    subroutine likeconsir(tau, lambda, n, ns, ni, tmin, tmax, alpha, beta, &
                           & spark, covmat, network, val)
    !log-likelihood for contact network model: SIR
    implicit none

    !Declarations
    integer, intent(in)           :: n, tmax, ns, ni, tmin
    integer, intent(in)           :: tau(n), lambda(n)          !inftime and infperiod
    double precision, intent(in)  :: alpha(ns), beta(ni), spark !parameters
    double precision, intent(in)  :: covmat(n,ns)               !covariates
    double precision, intent(in)  :: network(n,n,ni)            !contact network
    double precision, intent(out) :: val                        !result

    integer           :: i, j, t, k
    double precision  :: dx1, dx2, p1, p2, somega(n)

    somega = matmul(covmat,alpha) !susceptibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
        !infectious period
        if (tau(i)==(t+1)) then
          dx1 = 0.0d0
          do j = 1, n
          if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0) .AND. &
                            & (tau(j)+lambda(j) .GE. (t+1))) then
            do k = 1, ni
            dx1 = dx1 + beta(k) * network(i, j, k)
            end do
            end if
          end do
          p1 = 1.0d0 - exp(-(somega(i) * dx1 + spark))
          val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
          dx2 = 0.0d0
          do j = 1, n
            if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0) .AND.&
                       & (tau(j)+lambda(j) .GE. (t+1)))then
              do k = 1, ni
                dx2 = dx2 + beta(k) * network(i,j,k)
              end do
            end if
          end do
          p2 = exp(-(somega(i) * dx2 + spark))
          val = val + log(p2)  !result
        end if
      end do
    end do

    end subroutine likeconsir

    end module mcmcdata         !module ends



