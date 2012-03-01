! ============================================================================
! store_levels.f90
! ============================================================================
! Asigna valores iniciales a las variables necesarias para el RK4.

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


  subroutine store_levels

    use arrays
    
    implicit none
    
    if(grids_number.eq.1) then

       re_psi_p1 = re_psi1
       im_psi_p1 = im_psi1
       
       re_pi_p1 = re_pi1
       im_pi_p1 = im_pi1
       
       re_psi_a1 = re_psi_p1
       im_psi_a1 = im_psi_p1
       
       re_pi_a1 = re_pi_p1
       im_pi_a1 = im_pi_p1
       
    else if(grids_number.eq.2) then
       
       if(rank.lt.half_numtasks) then
          
          re_psi_p1 = re_psi1
          im_psi_p1 = im_psi1
          
          re_pi_p1 = re_pi1
          im_pi_p1 = im_pi1
          
          re_psi_a1 = re_psi_p1
          im_psi_a1 = im_psi_p1
          
          re_pi_a1 = re_pi_p1
          im_pi_a1 = im_pi_p1
          
       else  if(rank.ge.half_numtasks) then
          
          re_psi_p2 = re_psi2
          im_psi_p2 = im_psi2
          
          re_pi_p2 = re_pi2
          im_pi_p2 = im_pi2
          
          re_psi_a2 = re_psi_p2
          im_psi_a2 = im_psi_p2
          
          re_pi_a2 = re_pi_p2
          im_pi_a2 = im_pi_p2
       endif

    endif

  end subroutine store_levels
