! ============================================================================
! external_boundaries .f90
! ============================================================================
! Implementa las condiciones de frontera de la malla completa.

!     Copyright (C) 2012  Edison Montoya, eamonto@gmail.com

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Up to date: 29 Feb 2012					


  subroutine external_boundaries

    use arrays

    implicit none

    include 'mpif.h'

    real(double) :: aux1,aux2
    real(double) :: re_2der_psi,im_2der_psi,re_der_pi,im_der_pi


    if(boundary.eq.1) then  !free boundary

       if(grids_number.eq.1) then
       
          if(rank.eq.0) then

             aux1=1.0D5!0
             
             re_der_pi  = (-re_pi1(3)+four*re_pi1(2)-three*re_pi1(1))/(two*dx)
             im_der_pi  = (-im_pi1(3)+four*im_pi1(2)-three*im_pi1(1))/(two*dx)
          
             re_2der_psi  = (re_psi1(3)-two*re_psi1(2)+re_psi1(1))/(two*dx*dx)
             im_2der_psi  = (im_psi1(3)-two*im_psi1(2)+im_psi1(1))/(two*dx*dx)
             
             re_spi1(1) = aux1*(re_der_pi+re_2der_psi)/two
             im_spi1(1) = aux1*(im_der_pi+im_2der_psi)/two

          endif
          
          if(rank.eq.(numtasks-1)) then

             aux1=1.0D5!0
             
             re_der_pi  = (re_pi1(Nx-2)-four*re_pi1(Nx-1)+three*re_pi1(Nx))/(two*dx)
             im_der_pi  = (im_pi1(Nx-2)-four*im_pi1(Nx-1)+three*im_pi1(Nx))/(two*dx)
             
             re_2der_psi  = (re_psi1(Nx-2)-two*re_psi1(Nx-1)+re_psi1(Nx))/(two*dx*dx)
             im_2der_psi  = (im_psi1(Nx-2)-two*im_psi1(Nx-1)+im_psi1(Nx))/(two*dx*dx)
             
             re_spi1(Nx) = -aux1*(re_der_pi-re_2der_psi)/two
             im_spi1(Nx) = -aux1*(im_der_pi-im_2der_psi)/two

          endif
          

       else if(grids_number.eq.2) then
          
          if(rank.eq.0) then
             
             aux1=1.0D5!0
             
             re_der_pi  = (-re_pi1(3)+four*re_pi1(2)-three*re_pi1(1))/(two*dx)
             im_der_pi  = (-im_pi1(3)+four*im_pi1(2)-three*im_pi1(1))/(two*dx)
             
             re_2der_psi  = (re_psi1(3)-two*re_psi1(2)+re_psi1(1))/(two*dx*dx)
             im_2der_psi  = (im_psi1(3)-two*im_psi1(2)+im_psi1(1))/(two*dx*dx)
             
             re_spi1(1) = aux1*(re_der_pi+re_2der_psi)/two
             im_spi1(1) = aux1*(im_der_pi+im_2der_psi)/two
             
          else if(rank.eq.(half_numtasks-1)) then
             
             aux1=1.0D5!0
             
             re_der_pi  = (re_pi1(Nx-2)-four*re_pi1(Nx-1)+three*re_pi1(Nx))/(two*dx)
             im_der_pi  = (im_pi1(Nx-2)-four*im_pi1(Nx-1)+three*im_pi1(Nx))/(two*dx)
             
             re_2der_psi  = (re_psi1(Nx-2)-two*re_psi1(Nx-1)+re_psi1(Nx))/(two*dx*dx)
             im_2der_psi  = (im_psi1(Nx-2)-two*im_psi1(Nx-1)+im_psi1(Nx))/(two*dx*dx)
             
             re_spi1(Nx) = -aux1*(re_der_pi-re_2der_psi)/two
             im_spi1(Nx) = -aux1*(im_der_pi-im_2der_psi)/two
             
          else  if(rank.eq.half_numtasks) then
             
             aux1=1.0D5!0
             
             re_der_pi  = (-re_pi2(3)+four*re_pi2(2)-three*re_pi2(1))/(two*dx)
             im_der_pi  = (-im_pi2(3)+four*im_pi2(2)-three*im_pi2(1))/(two*dx)
             
             re_2der_psi  = (re_psi2(3)-two*re_psi2(2)+re_psi2(1))/(two*dx*dx)
             im_2der_psi  = (im_psi2(3)-two*im_psi2(2)+im_psi2(1))/(two*dx*dx)
             
             re_spi2(1) = aux1*(re_der_pi+re_2der_psi)/two
             im_spi2(1) = aux1*(im_der_pi+im_2der_psi)/two
             
          else  if(rank.eq.(numtasks-1)) then
             
             aux1=1.0D5!0
             
             re_der_pi  = (re_pi2(Nx-2)-four*re_pi2(Nx-1)+three*re_pi2(Nx))/(two*dx)
             im_der_pi  = (im_pi2(Nx-2)-four*im_pi2(Nx-1)+three*im_pi2(Nx))/(two*dx)
             
             re_2der_psi  = (re_psi2(Nx-2)-two*re_psi2(Nx-1)+re_psi2(Nx))/(two*dx*dx)
             im_2der_psi  = (im_psi2(Nx-2)-two*im_psi2(Nx-1)+im_psi2(Nx))/(two*dx*dx)
             
             re_spi2(Nx) = -aux1*(re_der_pi-re_2der_psi)/two
             im_spi2(Nx) = -aux1*(im_der_pi-im_2der_psi)/two
             
          endif

       endif

       
    else if(boundary.eq.2) then !boundaries from WDW

       if(grids_number.eq.1) then

          if(rank.eq.0) then

             aux1 = direction*sqrt(cons_pi*G/three)*(abs(x1(1)) -two*mu_0)/mu_0
             re_spi1(1) = aux1*( re_psi1(1)-re_psi1(2) )
             im_spi1(1) = aux1*( im_psi1(1)-im_psi1(2) )

          endif

          if(rank.eq.(numtasks-1)) then

             aux2 = direction*sqrt(cons_pi*G/three)*(abs(x1(Nx)) -two*mu_0)/mu_0
             re_spi1(Nx) = aux2*( re_psi1(Nx)-re_psi1(Nx-1) )
             im_spi1(Nx) = aux2*( im_psi1(Nx)-im_psi1(Nx-1) )
       
          endif


       else if(grids_number.eq.2) then
    
          if(rank.eq.0) then   

             aux1 = direction*sqrt(cons_pi*G/three)*(abs(x1(1)) -two*mu_0)/mu_0
             re_spi1(1) = aux1*( re_psi1(1)-re_psi1(2) )
             im_spi1(1) = aux1*( im_psi1(1)-im_psi1(2) )

          else if(rank.eq.(half_numtasks-1)) then

             aux2 = direction*sqrt(cons_pi*G/three)*(abs(x1(Nx)) -two*mu_0)/mu_0
             re_spi1(Nx) = aux2*( re_psi1(Nx)-re_psi1(Nx-1) )
             im_spi1(Nx) = aux2*( im_psi1(Nx)-im_psi1(Nx-1) )

          else  if(rank.eq.half_numtasks) then

             aux1 = direction*sqrt(cons_pi*G/three)*(abs(x2(1)) -two*mu_0)/mu_0
             re_spi2(1) = aux1*( re_psi2(1)-re_psi2(2) )
             im_spi2(1) = aux1*( im_psi2(1)-im_psi2(2) )
          
          else  if(rank.eq.(numtasks-1)) then

             aux2 = direction*sqrt(cons_pi*G/three)*(abs(x2(Nx)) -two*mu_0)/mu_0
             re_spi2(Nx) = aux2*( re_psi2(Nx)-re_psi2(Nx-1) )
             im_spi2(Nx) = aux2*( im_psi2(Nx)-im_psi2(Nx-1) )

          endif
       endif

    else
       print*,'Boundary condition not implemented!'
       stop
    endif

  end subroutine external_boundaries
