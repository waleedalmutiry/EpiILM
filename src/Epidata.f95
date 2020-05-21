!######################################################################
!# MODULE: subprograms
!#
!# AUTHORS:
!#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>,
!#         Waleed Almutiry <wkmtierie@qu.edu.sa>, and
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
!######################################################################

    module subprograms
    use, intrinsic :: iso_c_binding
    implicit none
    public :: dataxy, dataxysir, datacon, dataconsir

    contains

!######################################################################

    subroutine dataxy (x, y, n, tmin, tmax, ns, nt, ni, alpha, phi, beta, spark, &
                       & covmatsus, covmattrans, tau) bind(C, name="dataxy_")
    !Epidemic simulation under purely spatial model: SI

    external seedin
    external seedout
    external randomnumber

    integer (C_INT), intent(in)    :: n, tmax, ns, nt, ni, tmin
    real (C_DOUBLE), intent(in)    :: alpha(ns), phi(nt), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in)    :: covmatsus(n,ns), covmattrans(n,nt)               !covariates
    real (C_DOUBLE), intent(in)    :: x(n), y(n)                 !locations
    integer (C_INT), intent(inout) :: tau(n)                     !infection times

    integer (C_INT) :: i, j, t, A
    real (C_DOUBLE) :: u, dx, p
    real (C_DOUBLE) :: Somega(n), Tomega(n)


    !initialzing random seed
    call seedin()

    if ( ALL( tau == 0 ) )then
    call randomnumber(u)
    A = int(u * n) + 1  !Initialize infectious state
    tau(A) = tmin
    end if

    !Calculate the distance matrix
!    do i = 1, n
!        do j = i, n
!            eu(i,j) = sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
!            eu(j,i) = eu(i,j)
!        end do
!    end do

    Somega = matmul(covmatsus, alpha) !susceptibility function
    Tomega = matmul(covmattrans, phi) !transmissibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, (tmax - 1)
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0_c_double
              do j = 1,n
                if ((tau(j) .NE. 0) .and. (tau(j) .LE. t)) then
                  dx = dx + ((sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))**(-beta(ni)))*Tomega(j))
                end if
              end do
            p = 1.0_c_double - exp(-((Somega(i) * dx) + spark))
            call randomnumber(u)
            if (p .GT. u) then
              tau(i) = t + 1 !time at which individuals become infected
            end if
            end if
        end do
    end do

    call seedout()
    end subroutine dataxy

!######################################################################

    subroutine dataxysir (n, tmin, tmax, ns, nt, ni, alpha, phi, beta, spark, covmatsus, covmattrans, &
                    & lambda, x, y, tau, remt) bind(C, name="dataxysir_")
    !Epidemic simulation under purely spatial model: SIR

    external seedin
    external seedout
    external randomnumber

    integer (C_INT), intent(in)     :: n, tmax, ns, nt, ni, tmin
    integer (C_INT), intent(inout)  :: lambda(n)                  !infperiod
    real (C_DOUBLE), intent(in)     :: alpha(ns), phi(nt), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in)     :: covmatsus(n,ns), covmattrans(n,nt)               !covariates
    real (C_DOUBLE), intent(in)     :: x(n), y(n)                 !locations
    integer (C_INT), intent(inout)  :: remt(n), tau(n)            !removal time and inftime

    integer (C_INT)          :: i, j, t, A
    real (C_DOUBLE) :: u, dx, p
    real (C_DOUBLE) :: Somega(n), Tomega(n)


    !initialzing random seed
    call seedin()

    if ( ALL( tau .EQ. 0 ) )then
        call randomnumber(u)
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
!    do i = 1, n
!        do j = i, n
!            eu(i,j) = sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
!            eu(j,i) = eu(i,j)
!        end do
!    end do

    Somega = matmul(covmatsus, alpha) !susceptibility function
    Tomega = matmul(covmattrans, phi) !transmissibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, (tmax - 1)
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0_c_double
              do j = 1, n
                if ((tau(j) .NE. 0) ) then
                  if ((tau(j) .LE. t) .and. ((tau(j)+lambda(j)) .GT. t)) then
                    dx = dx + ((sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))**(-beta(ni)))*Tomega(j))
                  end if
                end if
              end do
            p = 1.0_c_double - exp(-((Somega(i) * dx) + spark))
            call randomnumber(u)
            if (p .GT. u) then
              tau(i) = t + 1  !time at which individuals become infected
              remt(i) = t + 1 + lambda(i)  !time at which individual is removed
            end if
            end if
        end do
    end do

    call seedout()
    end subroutine dataxysir

!######################################################################

    subroutine datacon(n, tmin, tmax, ns, nt, ni, alpha, phi, beta, &
                & spark, covmatsus, covmattrans, network, tau) bind(C, name="datacon_")
    !Epidemic simulation under contact network model: SI

    external seedin
    external seedout
    external randomnumber

    integer (C_INT), intent(in)    :: n, ns, nt, ni, tmax, tmin
    real (C_DOUBLE), intent(in)    :: alpha(ns), phi(nt), beta(ni), spark    !parameters
    real (C_DOUBLE), intent(in)    :: network(n,n,ni), covmatsus(n,ns), covmattrans(n,nt) !network and covariates
    integer (C_INT), intent(inout) :: tau(n)                        !inftime

    integer (C_INT) :: i, j, t, k, A
    real (C_DOUBLE) :: u, dx, p, Somega(n), Tomega(n)


    !initialzing random seed
    call seedin()

    if ( ALL( tau == 0 ) )then
        call randomnumber(u)
        A = int(u * n) + 1  !Initialize infectious state
        tau(A) = tmin
    end if

    Somega = matmul(covmatsus, alpha) !susceptibility function
    Tomega = matmul(covmattrans, phi) !transmissibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, (tmax - 1)
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0_c_double
              do j = 1, n
                if ((tau(j) .NE. 0) .and. (tau(j) .LE. t)) then
                  do k = 1, ni
                    dx = dx + ((beta(k) * network(i,j,k))*Tomega(j))
                  end do
                end if
              end do
            p = 1.0_c_double - exp(-((Somega(i) * dx) + spark))
            call randomnumber(u)
            if (p .GT. u) then
              tau(i) = t + 1  !time at which individuals become infected
            end if
            end if
        end do
    end do

    call seedout()
    end subroutine datacon

!######################################################################

    subroutine dataconsir(n, tmin, tmax, ns, nt, ni, lambda, alpha, phi, beta, &
                & spark, covmatsus, covmattrans, network, tau, remt) bind(C, name="dataconsir_")
    !Epidemic simulation under contact network model: SIR

    external seedin
    external seedout
    external randomnumber

    integer (C_INT), intent(in)     :: n, tmax, ns, nt, ni, tmin
    integer (C_INT), intent(in)     :: lambda(n)                     !infectious period
    real (C_DOUBLE), intent(in)     :: alpha(ns), phi(nt), beta(ni), spark    !parameters
    real (C_DOUBLE), intent(in)     :: network(n,n,ni), covmatsus(n,ns), covmattrans(n,nt) !network and covariates
    integer (C_INT), intent(inout)  :: tau(n)                        !infection times
    integer (C_INT), intent(inout)  :: remt(n)                       !removal time

    integer (C_INT) :: i, j, t, k, A
    real (C_DOUBLE) :: u, dx, p, Somega(n), Tomega(n)



    !initialzing random seed
    call seedin()

    if ( ALL( tau .EQ. 0 ) )then
        call randomnumber(u)
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

    Somega= matmul(covmatsus,alpha) !susceptibility function
    Tomega= matmul(covmattrans,phi) !transmissibility function

    !Calculate the probabilities of susceptible being exposed and update tau
    do t = tmin, (tmax - 1)
        do i = 1, n
            if (tau(i)==0) then
              dx = 0.0_c_double
              do j = 1, n
                if (tau(j) .NE. 0)  then
                  if ((tau(j) .LE. t).and. ((tau(j) + lambda(j)) .GT. t))then
                    do k = 1, ni
                      dx = dx +(( beta(k) * network(i,j,k))*Tomega(j))
                    end do
                  end if
                end if
              end do
            p = 1.0_c_double - exp(-((Somega(i) * dx) + spark))
            call randomnumber(u)
            if (p .GT. u) then
              tau(i) = t + 1 !time at which individuals become infected
              remt(i) = t + 1 + lambda(i) !time at which individual is removed
            end if
            end if
        end do
    end do

    call seedout()
    end subroutine dataconsir

    end module subprograms
