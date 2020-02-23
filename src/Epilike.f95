!######################################################################
!# MODULE: likesubprograms
!#
!# AUTHORS:
!#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca>,
!#         Waleed Almutiry <wkmtierie@qu.edu.sa>, and
!#         Rob Deardon <robert.deardon@ucalgary.ca>
!#
!# DESCRIPTION:
!#
!#     Calculates log-likelihood for spatial and network individual-level models with two
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
!#            like     .................. subroutine
!#            likesir  .................. subroutine
!#            likecon    .................. subroutine
!#            likeconsir .................. subroutine
!######################################################################

    module likesubprograms
    use ISO_C_BINDING
    implicit none
    public :: like, likesir, likecon, likeconsir

    contains

    subroutine like(x, y, tau, n, tmin, tmax, ns, nt, ni, alpha, phi, beta, &
                   & spark, covmatsus, covmattrans, val) bind(C, name="like_")
    !log-likelihood for purely spatial model: SI
    implicit none

    !Declarations
    integer (C_INT), intent(in) :: n, tau(n), tmax, ns, nt, ni, tmin
    real (C_DOUBLE), intent(in) :: alpha(ns), phi(nt), beta(ni), spark  !parameters
    real (C_DOUBLE), intent(in) :: x(n), y(n)                  !locations
    real (C_DOUBLE), intent(in) :: covmatsus(n, ns), covmattrans(n, nt)               !covariates
    real (C_DOUBLE), intent(out):: val                         !result

    integer (C_INT)          :: i, j, t
    real (C_DOUBLE) :: eu(n, n), Somega(n), Tomega(n)
    real (C_DOUBLE) :: dx1, dx2, p1, p2

    !Calculate the distance matrix
    do i = 1, n
      do j = i, n
        eu(i,j) = sqrt(((x(i)-x(j))**2) + ((y(i)-y(j))**2))
        eu(j,i) = eu(i,j)
      end do
    end do

    Somega = matmul(covmatsus,alpha) !susceptibility function
    Tomega = matmul(covmattrans,phi) !transmissibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
    !infectious period
        if (tau(i)==(t+1)) then
          dx1 = 0.0d0
          do j = 1, n
            if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0)) then
              dx1 = dx1 + ((eu(i,j)**(-beta(ni)))*Tomega(j))
            end if
          end do
        p1 = 1.0d0 - exp(-((Somega(i) * dx1) + spark))
        val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i) == 0)) then
          dx2 = 0.0d0
          do j = 1, n
            if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0))then
              dx2 = dx2 + ((eu(i,j)**(-beta(ni)))*Tomega(j))
            end if
          end do
        p2 = exp(-((Somega(i) * dx2) + spark))
        val = val + log(p2)  !result
        end if
      end do
    end do

    end subroutine like

!######################################################################

    subroutine likesir(x, y, tau, lambda, n, tmin, tmax, ns, nt, ni, alpha, &
                     & phi, beta, spark, covmatsus, covmattrans, val) bind(C, name="likesir_")
    !log-likelihood for purely spatial model: SIR
    implicit none

    !Declarations
    integer (C_INT), intent(in) :: n, tau(n), lambda(n), tmax, ns, nt, ni, tmin
    real (C_DOUBLE), intent(in) :: alpha(ns), beta(ni), phi(nt), spark  !parameters
    real (C_DOUBLE), intent(in) :: x(n), y(n)                  !locations
    real (C_DOUBLE), intent(in) :: covmatsus(n, ns), covmattrans(n, nt)               !covariates
    real (C_DOUBLE), intent(out):: val                         !result

    integer          :: i, j, t
    double precision :: eu(n,n), Somega(n), Tomega(n)
    double precision :: dx1, dx2, p1, p2

    !Calculate the distance matrix
    do i = 1, n
      do j = i, n
        eu(i,j) = sqrt(((x(i)-x(j))**2) + ((y(i)-y(j))**2))
        eu(j,i) = eu(i,j)
      end do
    end do

    Somega = matmul(covmatsus, alpha) !susceptibility function
    Tomega = matmul(covmattrans, phi) !transmissibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
    !infectious period
        if (tau(i)==(t+1)) then
          dx1 = 0.0d0
            do j = 1, n
              if (tau(j) .NE. 0) then
                if ((tau(j) .LT. (t+1)) .AND. (tau(j) + lambda(j) .GE. (t+1))) then
                  dx1 = dx1 + ((eu(i,j)**(-beta(ni)))*Tomega(j))
                end if
              end if
            end do
         p1 = 1.0d0 - exp(-((Somega(i) * dx1) + spark))
         val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
          dx2 = 0.0d0
            do j = 1, n
              if (tau(j) .NE. 0) then
                if ((tau(j) .LT. (t+1)) .AND. (tau(j)+lambda(j) .GE. (t+1)))then
                  dx2 = dx2 + ((eu(i,j)**(-beta(ni)))*Tomega(j))
                end if
              end if
            end do
          p2 = exp(-((Somega(i) * dx2) + spark))
          val = val + log(p2)  !result
        end if
      end do
    end do

    end subroutine likesir

!######################################################################

    subroutine likecon(tau, n, ns, nt, ni, tmin, tmax, alpha, phi, beta, &
                      & spark, covmatsus, covmattrans, network, val) bind(C, name="likecon_")
    !log-likelihood for contact network model: SI
    implicit none

    !Declarations
    integer (C_INT), intent(in)  :: n, tmax, ns, nt, ni, tmin
    integer (C_INT), intent(in)  :: tau(n)                     !infection times
    real (C_DOUBLE), intent(in)  :: alpha(ns), phi(nt), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in)  :: covmatsus(n,ns), covmattrans(n,nt)               !covariates
    real (C_DOUBLE), intent(in)  :: network(n,n,ni)            !contact network
    real (C_DOUBLE), intent(out) :: val                        !result

    integer          :: i, j, t, k
    double precision :: dx1, dx2, p1, p2, Somega(n), Tomega(n)

    Somega = matmul(covmatsus,alpha) !susceptibility function
    Tomega = matmul(covmattrans,phi) !transmissibility function

    val = 0.0d0
    do t = tmin, (tmax-1)
      do i = 1, n
      !infectious period
        if (tau(i)==(t+1)) then
         dx1 = 0.0d0
            do j = 1, n
              if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0)) then
                do k = 1, ni
                  dx1 = dx1 + ((beta(k) * network(i,j,k))*Tomega(j))
                end do
              end if
            end do
         p1 = 1.0d0 - exp(-((Somega(i) * dx1) + spark))
         val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
         dx2 = 0.0d0
            do j = 1, n
              if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0))then
                do k = 1, ni
                  dx2 = dx2 + ((beta(k) * network(i,j,k))*Tomega(j))
                end do
              end if
            end do
         p2 = exp(-((Somega(i) * dx2) + spark))
         val = val + log(p2)  !result
        end if
     end do
    end do

    end subroutine likecon

!######################################################################

    subroutine likeconsir(tau, lambda, n, ns, nt, ni, tmin, tmax, alpha, phi, beta, &
         & spark, covmatsus, covmattrans, network, val) bind(C, name="likeconsir_")
    !log-likelihood for contact network model: SIR
    implicit none

    !Declarations
    integer (C_INT), intent(in)  :: n, tmax, ns, nt, ni, tmin
    integer (C_INT), intent(in)  :: tau(n), lambda(n)          !inftime and infperiod
    real (C_DOUBLE), intent(in)  :: alpha(ns), phi(nt), beta(ni), spark !parameters
    real (C_DOUBLE), intent(in)  :: covmatsus(n,ns), covmattrans(n,nt)               !covariates
    real (C_DOUBLE), intent(in)  :: network(n,n,ni)            !contact network
    real (C_DOUBLE), intent(out) :: val                        !result

    integer           :: i, j, t, k
    double precision  :: dx1, dx2, p1, p2, Somega(n), Tomega(n)

    Somega = matmul(covmatsus,alpha) !susceptibility function
    Tomega = matmul(covmattrans,phi) !transmissibility function

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
                   dx1 = dx1 + ((beta(k) * network(i, j, k))*Tomega(j))
                 end do
              end if
            end do
         p1 = 1.0d0 - exp(-((Somega(i) * dx1) + spark))
         val = val + log(p1)
        end if
        !susceptible period
        if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
         dx2 = 0.0d0
            do j = 1, n
              if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0) .AND.&
                     & (tau(j)+lambda(j) .GE. (t+1)))then
                do k = 1, ni
                  dx2 = dx2 + ((beta(k) * network(i,j,k))*Tomega(j))
                end do
              end if
            end do
         p2 = exp(-((Somega(i) * dx2) + spark))
         val = val + log(p2)  !result
        end if
      end do
    end do

    end subroutine likeconsir

    end module likesubprograms
