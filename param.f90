! ===========================================================================
! param.f90
! ===========================================================================
! Global parameters for the physical system.

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


  module param

    implicit none
    
    INTEGER, PARAMETER :: double = SELECTED_REAL_KIND (13)
    
    !constants
    real(double),parameter :: mu_0 = 2.598076211D0!1.0D0 !
    real(double),parameter :: cons_pi = 3.141592654D0
    real(double),parameter :: G=1.0D0 !6.6732D-11!
    real(double),parameter :: hbar = 1.0D0!1.05459D-34!
    real(double),parameter :: c = 1.0D0
    real(double),parameter :: sgnp = -1.0D0
    real(double),parameter :: cons_K = 0.413602159600933D0
    real(double),parameter :: lambda = 2.584664094042158D0 !1.0D0
        
    !Dynamical Equations
    integer :: dyn_ec = 0 !0-> LQC !1-> sLQC_mu !2-> sLQC_v

    !output
    character(100) :: output_dir = "test_grid1"
    integer :: every_0D = 1    !at time every_0D*dt is write the time output
    integer :: every_1D = 10   !at time every_1D*dt is write the spacial output
    integer :: every_1D_grid = 1 !spacial step for output of grid points
    integer :: stdout = 1 !Output to Stardart out, 0--> OFF, 1--> ON 
    
    !grid
    integer,parameter :: grids_number = 2 
    integer :: Nx = 3000    !must be an even number for simetric grid
    real(double) :: epsilon =0.6*mu_0 ! 0.6*lambda !0.0D0 !

    !initial data
    integer :: method = 2 !method to find semi-classical state
    real(double) :: P_phi_star = sgnp*1000.0D0  !1000.0D0 !
    real(double) :: mu_star = 8000.0D0 !4000.0D0 !700.0D0 !
    real(double) :: phi_0 = 0.0D0
    real(double) :: delta_mu =300.0D0 !100.0D0 !
    
    !evolution
    integer :: boundary = 1   !1=wave 2=WDW
    real(double),parameter :: direction = -sgnp !'s' in WDW
    real(double) :: Total_time = direction*0.01D0 !direction*1.8D0!
    real(double) :: dt = 0.0001D0 !0.00001D0 !
    
    !numbers
    real(double) :: zero  = 0.0D0
    real(double) :: one   = 1.0D0
    real(double) :: two   = 2.0D0
    real(double) :: three = 3.0D0
    real(double) :: four  = 4.0D0
    
  end module param
