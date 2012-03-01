!============================================================================
!initialize.f90
!============================================================================
!Inicializa todos las variables a cero o valores predefinidos

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


  subroutine initialize

    use arrays
    
    implicit none

    include 'mpif.h'   

    !Observables
    norm = zero
    norm_sum = zero

    volume = zero
    volume_sum = zero

    re_P_phi = zero
    re_P_phi_sum = zero

    im_P_phi = zero
    im_P_phi_sum = zero

    !Time
    time = phi_0

    !Grids
    grid_center= 0 !(Nx)/two !<-- MAL CORREGIR

    !Auxiliars
    shift = 1.0D-15
    
    if (grids_number.eq.1)then

       x1 = zero
       if(dyn_ec.eq.1) sgn_x1 = zero
       
       !functions
       mod_psi1 = zero
       mod_pi1 = zero

       !Auxiliars
       auxiliar1 = zero
       auxiliar2 = zero

       B_mu1 = zero
       C_R1 = zero
       C_L1 = zero
       C_01 = zero
       
       re_psi1 = zero
       im_psi1 = zero
       re_pi1  = zero
       im_pi1 = zero
       
       re_psi_p1 = zero
       im_psi_p1 = zero
       re_pi_p1  = zero
       im_pi_p1  = zero
       
       re_psi_a1 = zero
       im_psi_a1 = zero
       re_pi_a1  = zero
       im_pi_a1  = zero
       
       re_spsi1 = zero
       im_spsi1 = zero
       re_spi1  = zero
       im_spi1  = zero

    else if(grids_number.eq.2) then

       if(rank.lt.half_numtasks) then
       
          x1 = zero
          if(dyn_ec.eq.1) sgn_x1 = zero
          
          !functions
          mod_psi1 = zero
          mod_pi1 = zero
          
          !Auxiliars
          auxiliar1 = zero
          auxiliar2 = zero
       
          B_mu1 = zero
          C_R1 = zero
          C_L1 = zero
          C_01 = zero
          
          re_psi1 = zero
          im_psi1 = zero
          re_pi1  = zero
          im_pi1 = zero
       
          re_psi_p1 = zero
          im_psi_p1 = zero
          re_pi_p1  = zero
          im_pi_p1  = zero
          
          re_psi_a1 = zero
          im_psi_a1 = zero
          re_pi_a1  = zero
          im_pi_a1  = zero
          
          re_spsi1 = zero
          im_spsi1 = zero
          re_spi1  = zero
          im_spi1  = zero
       
       
       else if(rank.ge.half_numtasks) then
          
          x2 = zero
          if(dyn_ec.eq.1) sgn_x2 = zero
          
          !functions
          mod_psi2 = zero
          mod_pi2 = zero
          
          !Auxiliars
          auxiliar1 = zero
          auxiliar2 = zero

          B_mu2 = zero
          C_R2 = zero
          C_L2 = zero
          C_02 = zero
          
          re_psi2 = zero
          im_psi2 = zero
          re_pi2  = zero
          im_pi2 = zero
          
          re_psi_p2 = zero
          im_psi_p2 = zero
          re_pi_p2  = zero
          im_pi_p2  = zero
          
          re_psi_a2 = zero
          im_psi_a2 = zero
          re_pi_a2  = zero
          im_pi_a2  = zero
          
          re_spsi2 = zero
          im_spsi2 = zero
          re_spi2  = zero
          im_spi2  = zero
          
       endif

    endif
 end subroutine initialize
