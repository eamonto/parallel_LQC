!============================================================================
!arrays.f90
!============================================================================
! Module with all the functions used.

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

  module arrays

    use param
    
    implicit none

    !Output
    character(100) :: log_file

    !Parallel
    integer :: rank,half_numtasks,numtasks,ierr  

    integer, allocatable, dimension (:) :: grid1_points_process
    integer, allocatable, dimension (:) :: grid2_points_process

    real(double), allocatable, dimension (:) :: grid1_length
    real(double), allocatable, dimension (:) :: grid2_length

    !Observables
    real(double) :: norm 
    real(double) :: norm_sum !For Parallel Communication

    real(double) :: volume
    real(double) :: volume_sum !For Parallel Communication

    real(double) :: re_P_phi
    real(double) :: re_P_phi_sum !For Parallel Communication

    real(double) :: im_P_phi
    real(double) :: im_P_phi_sum !For Parallel Communication

    !Time
    real(double) :: time

    !Grids
    real(double) :: length
    integer :: grid_center
    real(double) :: dx
    real(double), allocatable, dimension (:) :: x1
    real(double), allocatable, dimension (:) :: x2
    real(double), allocatable, dimension (:) :: sgn_x1
    real(double), allocatable, dimension (:) :: sgn_x2
    
    !Functions
    real(double), allocatable, dimension (:) :: mod_psi1
    real(double), allocatable, dimension (:) :: mod_psi2

    real(double), allocatable, dimension (:) :: mod_pi1
    real(double), allocatable, dimension (:) :: mod_pi2
    
    !Auxiliars
    real(double) :: shift

    real(double), allocatable, dimension (:) :: auxiliar1
    real(double), allocatable, dimension (:) :: auxiliar2

    real(double), allocatable, dimension (:) :: B_mu1
    real(double), allocatable, dimension (:) :: C_R1 !C+
    real(double), allocatable, dimension (:) :: C_L1 !C-
    real(double), allocatable, dimension (:) :: C_01 !C0

    real(double), allocatable, dimension (:) :: B_mu2
    real(double), allocatable, dimension (:) :: C_R2 !C+
    real(double), allocatable, dimension (:) :: C_L2 !C-
    real(double), allocatable, dimension (:) :: C_02 !C0
    
    real(double), allocatable, dimension (:) :: re_psi1,im_psi1
    real(double), allocatable, dimension (:) :: re_psi2,im_psi2

    real(double), allocatable, dimension (:) :: re_pi1, im_pi1
    real(double), allocatable, dimension (:) :: re_pi2, im_pi2
    
    real(double), allocatable, dimension (:) :: re_psi_p1,im_psi_p1
    real(double), allocatable, dimension (:) :: re_psi_p2,im_psi_p2

    real(double), allocatable, dimension (:) :: re_pi_p1, im_pi_p1
    real(double), allocatable, dimension (:) :: re_pi_p2, im_pi_p2
    
    real(double), allocatable, dimension (:) :: re_psi_a1,im_psi_a1
    real(double), allocatable, dimension (:) :: re_psi_a2,im_psi_a2

    real(double), allocatable, dimension (:) :: re_pi_a1, im_pi_a1
    real(double), allocatable, dimension (:) :: re_pi_a2, im_pi_a2
    
    real(double), allocatable, dimension (:) :: re_spsi1,im_spsi1
    real(double), allocatable, dimension (:) :: re_spsi2,im_spsi2

    real(double), allocatable, dimension (:) :: re_spi1, im_spi1
    real(double), allocatable, dimension (:) :: re_spi2, im_spi2

  end module arrays
