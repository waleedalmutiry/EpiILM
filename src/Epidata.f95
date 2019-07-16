!######################################################################
!# MODULE: subprograms
!#
!# AUTHORS:
!#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
!#         Rob Deardon <robert.deardon@ucalgary.ca>
!#
!# DESCRIPTION:
!#
!#     Simulates epidemics under spatial and network models with two
!#     disease types: SI and SIR
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
!#            dataxy     .................. subroutine
!#            dataxysir  .................. subroutine
!#            datacon    .................. subroutine
!#            dataconsir .................. subroutine
!#            initrandomseed  .............. subroutine
!######################################################################

    module subprograms
    use ISO_C_BINDING
    implicit none
    public :: dataxy, dataxysir, datacon, dataconsir

    contains

!######################################################################

    subroutine initrandomseed(tempseed)
    !The seed for the random number generation method
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

    subroutine dataxy (x, y, n, tmin, tmax, ns, ni, alpha, beta, spark, &
                       & covmat, tau, tempseed) bind(C, name="dataxy_")
    !Epidemic simulation under purely spatial model: SI
    implicit none

    integer (C_INT), intent(in)    :: n, tmax, ns, ni, tmin, tempseed
    real (C_DOUBLE), intent(in)    :: alpha(ns), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in)    :: covmat(n,ns)               !covariates
    real (C_DOUBLE), intent(in)    :: x(n), y(n)                 !locations
    integer (C_INT), intent(inout) :: tau(n)                     !infection times

    integer          :: i, j, t, A
    double precision :: u, dx, p
    double precision :: eu(n,n), Somega(n)


    !initialzing random seed
    if (tempseed .NE. 0) then
      call initrandomseed(tempseed)
    end if

    if ( ALL( tau == 0 ) )then
    call random_number(u)
    A = int(u * n) + 1  !Initialize infectious state
    tau(A) = tmin
    end if

    !Calculate the distance matrix
    do i = 1, n
        do j = i, n
            eu(i,j) = sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
            eu(j,i) = eu(i,j)
        end do
    end do

    Somega = matmul(covmat, alpha) !susceptibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, tmax
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0d0
              do j = 1,n
                if ((tau(j) .NE. 0) .and. (tau(j) .LE. t)) then
                  dx = dx + (eu(i,j)**(-beta(ni)))
                end if
              end do
            p = 1.0d0 - exp(-(Somega(i) * dx + spark))
            call random_number(u)
            if (p .GT. u) then
              tau(i) = t + 1 !time at which individuals become infected
            end if
            end if
        end do
    end do

    end subroutine dataxy

!######################################################################

    subroutine dataxysir (n, tmin, tmax, ns, ni, alpha, beta, spark, covmat, &
                    & lambda, x, y, tau, remt, tempseed) bind(C, name="dataxysir_")
    !Epidemic simulation under purely spatial model: SIR
    implicit none

    integer (C_INT), intent(in)     :: n, tmax, ns, ni, tmin, tempseed
    integer (C_INT), intent(inout)  :: lambda(n)                  !infperiod
    real (C_DOUBLE), intent(in)     :: alpha(ns), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in)     :: covmat(n,ns)               !covariates
    real (C_DOUBLE), intent(in)     :: x(n), y(n)                 !locations
    integer (C_INT), intent(inout)  :: remt(n), tau(n)            !removal time and inftime

    integer          :: i, j, t, A
    double precision :: u, dx, p
    double precision :: eu(n,n), Somega(n)


    !initialzing random seed
    if (tempseed .NE. 0) then
      call initrandomseed(tempseed)
    end if

    if ( ALL( tau .EQ. 0 ) )then
        call random_number(u)
        A = int(u * n) + 1  !Initialize infectious state
        tau(A) = tmin
        remt(A) = tau(A) + lambda(A)
    else
        do i = 1, n
            if (tau(i) .NE. 0) then
                remt(i) = tau(i) + lambda(i)
            else
                 remt(i) = tau(i)
            end if
        end do
   end if

    !Calculate the distance matrix
    do i = 1, n
        do j = i, n
            eu(i,j) = sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
            eu(j,i) = eu(i,j)
        end do
    end do

    Somega = matmul(covmat, alpha) !susceptibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, tmax
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0d0
              do j = 1, n
                if ((tau(j) .NE. 0) ) then
                  if ((tau(j) .LE. t) .and. ((tau(j)+lambda(j)) .GT. t)) then
                    dx = dx + (eu(i,j)**(-beta(ni)))
                  end if
                end if
              end do
            p = 1.0d0 - exp(-(Somega(i) * dx + spark))
            call random_number(u)
            if (p .GT. u) then
              tau(i) = t + 1  !time at which individuals become infected
              remt(i) = t + 1 + lambda(i)  !time at which individual is removed
            end if
            end if
        end do
    end do

    end subroutine dataxysir

!######################################################################

    subroutine datacon(n, tmin, tmax, ns, ni, alpha, beta, &
                & spark, covmat, network, tau, tempseed) bind(C, name="datacon_")
    !Epidemic simulation under contact network model: SI
    implicit none

    integer (C_INT), intent(in)    :: n, ns, ni, tmax, tmin, tempseed
    real (C_DOUBLE), intent(in)    :: alpha(ns), beta(ni), spark    !parameters
    real (C_DOUBLE), intent(in)    :: network(n,n,ni), covmat(n,ns) !network and covariates
    integer (C_INT), intent(inout) :: tau(n)                        !inftime

    integer          :: i, j, t, k, A
    double precision :: u, dx, p, Somega(n)


    !initialzing random seed
    if (tempseed .NE. 0) then
      call initrandomseed(tempseed)
    end if

    if ( ALL( tau == 0 ) )then
        call random_number(u)
        A = int(u * n) + 1  !Initialize infectious state
        tau(A) = tmin
    end if

    Somega = matmul(covmat, alpha) !susceptibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, tmax
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0d0
              do j = 1, n
                if ((tau(j) .NE. 0) .and. (tau(j) .LE. t)) then
                  do k = 1, ni
                    dx = dx + beta(k) * network(i,j,k)
                  end do
                end if
              end do
            p = 1.0d0 - exp(-(Somega(i) * dx + spark))
            call random_number(u)
            if (p .GT. u) then
              tau(i) = t + 1  !time at which individuals become infected
            end if
            end if
        end do
    end do

    end subroutine datacon

!######################################################################

    subroutine dataconsir(n, tmin, tmax, ns, ni, lambda, alpha, beta, &
                & spark, covmat, network, tau, remt, tempseed) bind(C, name="dataconsir_")
    !Epidemic simulation under contact network model: SIR
    implicit none

    integer (C_INT), intent(in)     :: n, tmax, ns, ni, tmin, tempseed
    integer (C_INT), intent(in)     :: lambda(n)                     !infectious period
    real (C_DOUBLE), intent(in)     :: alpha(ns), beta(ni), spark    !parameters
    real (C_DOUBLE), intent(in)     :: network(n,n,ni), covmat(n,ns) !network and covariates
    integer (C_INT), intent(inout)  :: tau(n)                        !infection times
    integer (C_INT), intent(inout)  :: remt(n)                       !removal time

    integer          :: i, j, t, k, A
    double precision :: u, dx, p, Somega(n)



    !initialzing random seed
    if (tempseed .NE. 0) then
      call initrandomseed(tempseed)
    end if

    if ( ALL( tau .EQ. 0 ) )then
        call random_number(u)
        A = int(u * n) + 1  !Initialize infectious state
        tau(A) = tmin
        remt(A) = tau(A) + lambda(A)
    else
        do i = 1, n
            if (tau(i) .NE. 0) then
                remt(i) = tau(i) + lambda(i)
            else
                remt(i) = tau(i)
            end if
        end do
    end if

    Somega= matmul(covmat,alpha) !susceptibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, tmax
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0d0
              do j = 1, n
                if (tau(j) .NE. 0)  then
                  if ((tau(j) .LE. t).and. ((tau(j) + lambda(j)) .GT. t))then
                    do k = 1, ni
                      dx = dx + beta(k) * network(i,j,k)
                    end do
                  end if
                end if
              end do
            p = 1.0d0 - exp(-(Somega(i) * dx + spark))
            call random_number(u)
            if (p .GT. u) then
              tau(i) = t + 1 !time at which individuals become infected
              remt(i) = t + 1 + lambda(i) !time at which individual is removed
            end if
            end if
        end do
    end do

    end subroutine dataconsir

    end module subprograms







