! ============================================================================
! observables.f90
! ============================================================================
! Calcula el valor de los oservables: volume, P_phi (y la norma).
! Cada proceso calcula su parte y luego se unifican los valores en el
! proceso cero.

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


subroutine observables
  
  use arrays
  
  implicit none
  
  include 'mpif.h'
  
  integer :: i
  
  norm = zero
  norm_sum = zero
  
  volume = zero
  volume_sum = zero
  
  re_P_phi = zero
  re_P_phi_sum = zero
  
  im_P_phi = zero
  im_P_phi_sum = zero
  
  !Queda faltando la contribucion de los extremos de la malla a los observables
  
  if(dyn_ec.eq.0) then
     
     if(grids_number.eq.1) then
        
        mod_psi1 = re_psi1**two + im_psi1**two
        mod_pi1  = re_pi1**two  + im_pi1**two
        
        do i=2,Nx-1
           norm = norm + B_mu1(i)*mod_psi1(i)*dx
           
           volume = volume + B_mu1(i)*abs(x1(i))*mod_psi1(i)*dx
           
           re_P_phi = re_P_phi + B_mu1(i)*(re_psi1(i)*im_pi1(i) & 
                -im_psi1(i)*re_pi1(i))*dx
           
           im_P_phi = im_P_phi + B_mu1(i)*(-re_psi1(i)*re_pi1(i) & 
                -im_psi1(i)*im_pi1(i))*dx
        enddo
        
        call MPI_REDUCE(norm,norm_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(volume,volume_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(re_P_phi,re_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(im_P_phi,im_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(rank.eq.0) then
           volume_sum=volume_sum/norm_sum
           re_P_phi_sum = re_P_phi_sum/norm_sum
           im_P_phi_sum = im_P_phi_sum/norm_sum
        endif
          
     
     else if(grids_number.eq.2) then
          
        if(rank.lt.half_numtasks) then
           
           mod_psi1 = re_psi1**two + im_psi1**two
           mod_pi1  = re_pi1**two  + im_pi1**two
           
           do i=2,Nx-1
              norm = norm + B_mu1(i)*mod_psi1(i)*dx
              
              volume = volume + B_mu1(i)*abs(x1(i))*mod_psi1(i)*dx
              
              re_P_phi = re_P_phi + B_mu1(i)*(re_psi1(i)*im_pi1(i) & 
                   -im_psi1(i)*re_pi1(i))*dx
              
              im_P_phi = im_P_phi + B_mu1(i)*(-re_psi1(i)*re_pi1(i) & 
                   -im_psi1(i)*im_pi1(i))*dx
           enddo
           
        else if(rank.ge.half_numtasks) then
           
           mod_psi2 = re_psi2**two + im_psi2**two
           mod_pi2  = re_pi2**two  + im_pi2**two
           
           do i=2,Nx-1
              norm = norm + B_mu2(i)*mod_psi2(i)*dx
              
              volume = volume + B_mu2(i)*abs(x2(i))*mod_psi2(i)*dx
              
              re_P_phi = re_P_phi + B_mu2(i)*(re_psi2(i)*im_pi2(i) & 
                   -im_psi2(i)*re_pi2(i))*dx
              
              im_P_phi = im_P_phi + B_mu2(i)*(-re_psi2(i)*re_pi2(i) & 
                   -im_psi2(i)*im_pi2(i))*dx
           enddo
             
        endif

        call MPI_REDUCE(norm,norm_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(volume,volume_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(re_P_phi,re_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(im_P_phi,im_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(rank.eq.0) then
           volume_sum=volume_sum/norm_sum
           re_P_phi_sum = re_P_phi_sum/norm_sum
           im_P_phi_sum = im_P_phi_sum/norm_sum
        endif
          
     end if
       
  else if(dyn_ec.eq.1) then
    
     !NO IMPLEMENTADO TODAVIA
     
  else if(dyn_ec.eq.2) then
       
     if(grids_number.eq.1) then

        mod_psi1 = re_psi1**two + im_psi1**two
        mod_pi1  = re_pi1**two  + im_pi1**two
          
        !For no make a zero division
        if(epsilon.eq.zero) B_mu1(grid_center)=one
        
        do i=2,Nx-1
           norm = norm + mod_psi1(i)*dx/abs(B_mu1(i))
           
           volume = volume + mod_psi1(i)*dx
           
           re_P_phi = re_P_phi + (re_psi1(i)*im_pi1(i) & 
                -im_psi1(i)*re_pi1(i))*dx/abs(B_mu1(i))
           
           im_P_phi = im_P_phi + (-re_psi1(i)*re_pi1(i) & 
                -im_psi1(i)*im_pi1(i))*dx/abs(B_mu1(i))
        enddo
        
        norm=norm*lambda/cons_pi
        volume=volume*(lambda**three)/(cons_pi*sqrt(three))
        re_P_phi = re_P_phi*hbar*lambda/cons_pi
        im_P_phi = im_P_phi*hbar*lambda/cons_pi

        call MPI_REDUCE(norm,norm_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(volume,volume_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(re_P_phi,re_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(im_P_phi,im_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        
        if(rank.eq.0) then
           volume_sum=volume_sum/norm_sum
           re_P_phi_sum = re_P_phi_sum/norm_sum
           im_P_phi_sum = im_P_phi_sum/norm_sum
        endif

        ! Restauration of original value
        if(epsilon.eq.zero) B_mu1(grid_center)=zero
          
     else if(grids_number.eq.2) then
          
        if(rank.lt.half_numtasks) then
           
           mod_psi1 = re_psi1**two + im_psi1**two
           mod_pi1  = re_pi1**two  + im_pi1**two
           
           !Creo que no es necesario, no hay dos mallas con epsilon=cero
           !For no make a zero division
           if(epsilon.eq.zero) B_mu1(grid_center)=one
           
           do i=2,Nx-1
              norm = norm + mod_psi1(i)*dx/abs(B_mu1(i))
              
              volume = volume + mod_psi1(i)*dx
                
              re_P_phi = re_P_phi + (re_psi1(i)*im_pi1(i) & 
                   -im_psi1(i)*re_pi1(i))*dx/abs(B_mu1(i))
              
              im_P_phi = im_P_phi + (-re_psi1(i)*re_pi1(i) & 
                   -im_psi1(i)*im_pi1(i))*dx/abs(B_mu1(i))
           enddo
           
        else if(rank.ge.half_numtasks) then
             
           mod_psi2 = re_psi2**two + im_psi2**two
           mod_pi2  = re_pi2**two  + im_pi2**two
           
           do i=2,Nx-1
              norm = norm + mod_psi2(i)*dx/abs(B_mu2(i))
              
              volume = volume + mod_psi2(i)*dx
              
              re_P_phi = re_P_phi + (re_psi2(i)*im_pi2(i) & 
                   -im_psi2(i)*re_pi2(i))*dx/abs(B_mu2(i))
              
              im_P_phi = im_P_phi + (-re_psi2(i)*re_pi2(i) & 
                   -im_psi2(i)*im_pi2(i))*dx/abs(B_mu2(i))
           enddo
           
        endif
        
        norm=norm*lambda/cons_pi
        volume=volume*(lambda**three)/(cons_pi*sqrt(three))
        re_P_phi = re_P_phi*hbar*lambda/cons_pi
        im_P_phi = im_P_phi*hbar*lambda/cons_pi
        
        call MPI_REDUCE(norm,norm_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(volume,volume_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(re_P_phi,re_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(im_P_phi,im_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       
        if(rank.eq.0) then
           volume_sum=volume_sum/norm_sum
           re_P_phi_sum = re_P_phi_sum/norm_sum
           im_P_phi_sum = im_P_phi_sum/norm_sum
        endif

        !Creo que no es necesario, no hay dos mallas con epsilon=cero
        ! Restauration of original value
        if(epsilon.eq.zero) B_mu1(grid_center)=zero

     endif
    
  endif
    
end subroutine observables

!!$       call MPI_ALLREDUCE(norm,norm_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$ !      call MPI_ALLREDUCE(volume,volume_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$       call MPI_ALLREDUCE(re_P_phi,re_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$       call MPI_ALLREDUCE(im_P_phi,im_P_phi_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$       
!!$       if(numtasks.gt.1) then
!!$          if(rank.eq.0) then
!!$             
!!$!             print *,'volume procesador 0',volume
!!$             volume_sum=volume
!!$             
!!$             do i=1,numtasks-1
!!$                call MPI_RECV(volume,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,stat,ierr) 
!!$ !               print *,'volume procesador',i,volume
!!$                volume_sum=volume_sum+volume
!!$             enddo
!!$  !           print *,'suma volume',volume_sum
!!$             volume_sum=volume_sum/norm_sum
!!$   !          print *,'suma volume real',volume_sum
!!$          else
!!$             call MPI_SEND(volume,1,MPI_DOUBLE_PRECISION,0,rank,MPI_COMM_WORLD,ierr)
!!$          endif
!!$       endif
       



