!######################################################################
!# MODULE: subprogramr0
!#
!# AUTHORS:
!#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>,
!#         Waleed Almutiry <wkmtierie@qu.edu.sa>, and
!#         Rob Deardon <robert.deardon@ucalgary.ca>
!#
!# DESCRIPTION:
!#
!#     Calculates the basic reproduction number for the specified SIR model and data
!#
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
!#            rxysir   ..................... subroutine
!#            rconsir  ..................... subroutine
!######################################################################

    module  subprogramr0
    use ISO_C_BINDING
    implicit none
    public :: rxysir,rconsir

    contains

!######################################################################

    subroutine rxysir (n, tmax, ns, nt, ni, suspar, transpar, beta, spark, &
                      & covmatsus, covmattrans, lambda, x, y, sim, &
                      & val, countinf)  bind(C, name="rxysir_")
    !Calculates the basic reproduction number under the spatial model
    implicit none

    external seedin
    external seedout
    external randomnumber

    integer (C_INT), intent(in)  :: n, tmax, ns, ni, nt, sim
    integer (C_INT), intent(in) :: lambda(n)                  !infperiod

    real (C_DOUBLE), intent(in) :: suspar(ns), transpar(nt), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in) :: covmatsus(n,ns), covmattrans(n,ns)        !susceptibility covariate
    real (C_DOUBLE), intent(in) :: x(n), y(n)                 !locations

    real (C_DOUBLE), intent(inout) :: val                     !result
    integer (C_INT), intent(out)   :: countinf(sim)           !result

    integer          :: i, j, t, tau(n), A
    double precision :: u, dx, p
    double precision :: eu(n,n), Somega(n), Tomega(n)

    !initialzing random seed
    call seedin()
 

    !Calculate the distance matrix
    do i = 1,n
        do j = i,n
            eu(i,j) = sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
            eu(j,i) = eu(i,j)
        end do
    end do

    Somega = matmul(covmatsus, suspar) !susceptibility function
    Tomega = matmul(covmattrans, transpar) !transmissibility function

    do j = 1,sim  !simulation starts here

    call randomnumber(u)
    A = int(u * n) + 1  !random initialization for first infection
    do i = 1,n
        tau(i) = 0
    end do
    tau(A) = 1

    !Calculate the initial probability of susceptible being exposed and update tau
    do t = 1,tmax
        do i = 1,n
          dx = 0.0d0
          if (tau(i)==0) then
              if ((tau(A) .LE. t) .and. ((tau(A)+lambda(A)) .GT. t)) then
                dx = (eu(i,A)**(-beta(ni)))*Tomega(A)
                p = 1.0d0 - exp(-((Somega(i)*dx)+spark))
                call randomnumber(u)
                if (p .GT. u) then
                  tau(i) = t + 1  !time at which individual become infected/infectious
                end if
              end if
          end if
        end do
    end do

    countinf(j) = 0
    do i = 1,n
    if (tau(i) .NE. 0) then
        countinf(j) = countinf(j) + 1
    end if
    end do

    countinf(j) = countinf(j)-1

    end do ! simulation ends

    val = sum(countinf)/dble(sim) !result

    call seedout()
    end subroutine rxysir

!######################################################################

    subroutine rconsir(n, tmax, ns, nt, ni, lambda, suspar, transpar, beta, &
                       & spark, covmatsus, covmattrans,network, sim, &
                       & val, countinf) bind(C, name="rconsir_")
    !Calculates the basic reproduction number under the contact network model
    implicit none

    external seedin
    external seedout
    external randomnumber

    integer (C_INT),intent(in)  :: n, tmax, ns, nt, ni, sim
    integer (C_INT),intent(in)  :: lambda(n)                   !infperiod
    real (C_DOUBLE), intent(in) :: suspar(ns), transpar(nt), beta(ni), spark  !parameter
    real (C_DOUBLE), intent(in) :: network(n, n, ni)           !contact network
    real (C_DOUBLE), intent(in) :: covmatsus(n, ns), covmattrans(n, nt)      !susceptibility covariates

    real (C_DOUBLE), intent(inout) :: val                      !result
    integer (C_INT),intent(out)    :: countinf(sim)            !result

    integer          :: i, j, t, k, tau(n), A
    double precision :: u, dx, p, Somega(n), Tomega(n)

    !initialzing random seed
    call seedin()
    

    Somega = matmul(covmatsus,suspar) !susceptibility function
    Tomega = matmul(covmattrans,transpar) !transmissibility function

    do j = 1, sim  !simulation starts

    call randomnumber(u)
    A = int(u * n) + 1 !random initialization for first infection
    do i = 1, n
        tau(i) = 0
    end do
    tau(A) = 1
    !Calculate the initial probability of susceptible being exposed and update tau
    do t = 1, tmax
        do i = 1, n
          dx = 0.0d0
          if (tau(i)==0) then
             if ((tau(A) .LE. t) .and. ((tau(A) + lambda(A)) .GT. t)) then
                do k = 1, ni
                  dx = dx + ((beta(k) * network(i, A, k))*Tomega(A))
                end do
                p = 1.0d0 - (exp(-((Somega(i) * dx) + spark)))
                call randomnumber(u)
                if (p .GT. u) then
                  tau(i) = t + 1
                end if
             end if
          end if
        end do
    end do

    countinf(j) = 0
    do i = 1, n
    if (tau(i) .NE. 0) then
      countinf(j) = countinf(j) + 1
    end if
    end do

    countinf(j) = countinf(j) - 1

    end do ! simulation ends

    val = sum(countinf)/dble(sim) !result

    call seedout()
    end subroutine rconsir

    end module subprogramr0
