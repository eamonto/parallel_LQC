! ============================================================================
! evolution.f90
! ============================================================================
! Coordina la evolución de las variables que se estan integrando en cada
! paso del RK4.

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


  subroutine evolution(k)

    use arrays

    implicit none

    include 'mpif.h'

    integer :: k
    real(double) :: dtw,weight

    if (k.eq.1) then

       dtw = dt/two
       weight = dt/6.0D0
       call step_rk4(dtw,weight)

    else if (k.eq.2) then

       dtw = dt/two
       weight = dt/three
       call step_rk4(dtw,weight)

    else if (k.eq.3) then

       dtw = dt
       weight = dt/three
       call step_rk4(dtw,weight)

    else if (k.eq.4) then

       weight = dt/6.0D0

       if(grids_number.eq.1) then
          re_psi1 = re_psi_a1 + weight*re_spsi1
          im_psi1 = im_psi_a1 + weight*im_spsi1
          
          re_pi1 = re_pi_a1 + weight*re_spi1
          im_pi1 = im_pi_a1 + weight*im_spi1
          
       else if(grids_number.eq.2) then
          
          if(rank.lt.half_numtasks) then
             
             re_psi1 = re_psi_a1 + weight*re_spsi1
             im_psi1 = im_psi_a1 + weight*im_spsi1
             
             re_pi1 = re_pi_a1 + weight*re_spi1
             im_pi1 = im_pi_a1 + weight*im_spi1
             
          else  if(rank.ge.half_numtasks) then
             
             re_psi2 = re_psi_a2 + weight*re_spsi2
             im_psi2 = im_psi_a2 + weight*im_spsi2
             
             re_pi2 = re_pi_a2 + weight*re_spi2
             im_pi2 = im_pi_a2 + weight*im_spi2
          endif
       endif
    end if

  end subroutine evolution
