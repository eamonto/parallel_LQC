! ============================================================================
! sources.f90
! ============================================================================
! Asigna los valores (dados por las ecuaciones de LQC) a las variables 
! que se van a integrar.

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


  subroutine sources

    use arrays

    implicit none

    integer :: i
    real(double) :: aux
    
    if(dyn_ec.eq.0) then

       if(grids_number.eq.1) then
          
          re_spsi1 = re_pi1
          im_spsi1 = im_pi1
          
          do i=2,Nx-1
             re_spi1(i) = (one/B_mu1(i))*(C_R1(i)*re_psi1(i+1)+C_01(i)*re_psi1(i)+C_L1(i)*re_psi1(i-1))
             im_spi1(i) = (one/B_mu1(i))*(C_R1(i)*im_psi1(i+1)+C_01(i)*im_psi1(i)+C_L1(i)*im_psi1(i-1))
          enddo
          

          !***ESTO HAY QUE CAMBIARLO***
          !THIS IS BECAUSE THERE ARE A ZERO DIVISION
          if(epsilon.eq.zero) then
             re_spi1(grid_center)=shift
             im_spi1(grid_center)=shift
          endif
          !***ESTO HAY QUE CAMBIARLO***

       
       else if(grids_number.eq.2) then

          if(rank.lt.half_numtasks) then
             
             re_spsi1 = re_pi1
             im_spsi1 = im_pi1
             
             do i=2,Nx-1
                re_spi1(i) = (one/B_mu1(i))*(C_R1(i)*re_psi1(i+1)+C_01(i)*re_psi1(i)+C_L1(i)*re_psi1(i-1))
                im_spi1(i) = (one/B_mu1(i))*(C_R1(i)*im_psi1(i+1)+C_01(i)*im_psi1(i)+C_L1(i)*im_psi1(i-1))
             enddo

             !***ESTO NO ES NECESARIO, NO HAY DOS MALLAS CON EPSILON=CERO***
             !THIS IS BECAUSE THERE ARE A ZERO DIVISION
             if(epsilon.eq.zero) then
                re_spi1(grid_center)=shift
                im_spi1(grid_center)=shift
             endif
             !***ESTO NO ES NECESARIO, NO HAY DOS MALLAS CON EPSILON=CERO***

          else  if(rank.ge.half_numtasks) then
          
             re_spsi2 = re_pi2
             im_spsi2 = im_pi2
          
             do i=2,Nx-1
                re_spi2(i) = (one/B_mu2(i))*(C_R2(i)*re_psi2(i+1)+C_02(i)*re_psi2(i)+C_L2(i)*re_psi2(i-1))
                im_spi2(i) = (one/B_mu2(i))*(C_R2(i)*im_psi2(i+1)+C_02(i)*im_psi2(i)+C_L2(i)*im_psi2(i-1))
             enddo
          
          endif
    
       endif
   
    else if(dyn_ec.eq.1) then

       !FALTA POR IMPLEMENTAR

    else if(dyn_ec.eq.2) then

       if(grids_number.eq.1) then

          aux=three*cons_pi*G/(four*lambda*lambda)

          re_spsi1 = re_pi1
          im_spsi1 = im_pi1

          do i=2,Nx-1
             re_spi1(i) = aux*B_mu1(i)*(C_R1(i)*re_psi1(i+1)+C_01(i)*re_psi1(i)+C_L1(i)*re_psi1(i-1))
             im_spi1(i) = aux*B_mu1(i)*(C_R1(i)*im_psi1(i+1)+C_01(i)*im_psi1(i)+C_L1(i)*im_psi1(i-1))
          enddo
       
          !THIS IS BECAUSE 'shift=1.0D-15' IS THE NUMERICAL ZERO 
          if(epsilon.eq.zero) then
             re_spi1(grid_center)=shift
             im_spi1(grid_center)=shift
          endif

       else if(grids_number.eq.2) then
          
          if(rank.lt.half_numtasks) then
             
             aux=three*cons_pi*G/(four*lambda*lambda)
             
             re_spsi1 = re_pi1
             im_spsi1 = im_pi1
             
             do i=2,Nx-1
                re_spi1(i) = aux*B_mu1(i)*(C_R1(i)*re_psi1(i+1)+C_01(i)*re_psi1(i)+C_L1(i)*re_psi1(i-1))
                im_spi1(i) = aux*B_mu1(i)*(C_R1(i)*im_psi1(i+1)+C_01(i)*im_psi1(i)+C_L1(i)*im_psi1(i-1))
             enddo
             
             !***ESTO NO ES NECESARIO, NO HAY DOS MALLAS CON EPSILON=CERO***
             !THIS IS BECAUSE 'shift=1.0D-15' IS THE NUMERICAL ZERO 
             if(epsilon.eq.zero) then
                re_spi1(grid_center)=shift
                im_spi1(grid_center)=shift
             endif
             !***ESTO NO ES NECESARIO, NO HAY DOS MALLAS CON EPSILON=CERO***

          else  if(rank.ge.half_numtasks) then

             aux=three*cons_pi*G/(four*lambda*lambda)
             
             re_spsi2 = re_pi2
             im_spsi2 = im_pi2
             
             do i=2,Nx-1
                re_spi2(i) = aux*B_mu2(i)*(C_R2(i)*re_psi2(i+1)+C_02(i)*re_psi2(i)+C_L2(i)*re_psi2(i-1))
                im_spi2(i) = aux*B_mu2(i)*(C_R2(i)*im_psi2(i+1)+C_02(i)*im_psi2(i)+C_L2(i)*im_psi2(i-1))
             enddo

          endif

       endif

    endif

  end subroutine sources
