!============================================================================
!allocate.f90
!============================================================================
!Alocación de la memoria para cada procesador

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


  subroutine allocate

    use arrays
    
    implicit none

    include 'mpif.h'

    integer :: i

    if (grids_number.eq.1)then
     
       !Numero de puntos a utilizar en los procesadores
       if(rank.eq.0) then
          Nx=grid1_points_process(rank)+2
       else if(rank.eq.(numtasks-1)) then
          Nx=grid1_points_process(rank)+1
       else
          Nx=grid1_points_process(rank)+2
       endif
       !Numero de puntos a utilizar en los procesadores
       
       if (numtasks.eq.1) Nx=grid1_points_process(0)+1

       !Output to the log file           
       if(rank.eq.0) then
          open(200,file='log_'//trim(output_dir),form='formatted',status='old',position='append')          
          write(200,"()")
          write(200,"(I8A45I8)") Nx,'puntos reales de la malla 1 para el proceso',rank
          
          do i=1,numtasks-2
             write(200,"(I8A45I8)") grid1_points_process(i)+2,'puntos reales de la malla 1 para el proceso',i
          enddo
          
          write(200,"(I8A45I8)") grid1_points_process(numtasks-1)+1,'puntos reales de la malla 1 para el proceso',numtasks-1
          write(200,"()")
          close(200)
       endif
       !Output to the log file           


       !Allocation of memory for first grid
       allocate(x1(1:Nx))
       if(dyn_ec.eq.1) allocate(sgn_x1(1:Nx))
       
       allocate(mod_psi1(1:Nx))
       allocate(mod_pi1(1:Nx))
       
       allocate(auxiliar1(1:Nx))
       allocate(auxiliar2(1:Nx))!se necesita para el metodo 1

       allocate(B_mu1(1:Nx))
       allocate(C_R1(1:Nx))
       allocate(C_L1(1:Nx))
       allocate(C_01(1:Nx))
       
       allocate(re_psi1(1:Nx))
       allocate(im_psi1(1:Nx))
       
       allocate(re_pi1(1:Nx))
       allocate(im_pi1(1:Nx))

       allocate(re_psi_p1(1:Nx))
       allocate(im_psi_p1(1:Nx))
       allocate(re_pi_p1(1:Nx))
       allocate(im_pi_p1(1:Nx))
       
       allocate(re_psi_a1(1:Nx))
       allocate(im_psi_a1(1:Nx))
       allocate(re_pi_a1(1:Nx))
       allocate(im_pi_a1(1:Nx))
    
       allocate(re_spsi1(1:Nx))
       allocate(im_spsi1(1:Nx))
       allocate(re_spi1(1:Nx))
       allocate(im_spi1(1:Nx))


    else if(grids_number.eq.2) then

       if(rank.lt.half_numtasks) then
          
          !Numero de puntos a utilizar en los procesadores
          if(rank.eq.0) then
             Nx=grid1_points_process(rank)+2
             
          else if(rank.eq.(half_numtasks-1)) then
             Nx=grid1_points_process(rank)+1
             
          else
             Nx=grid1_points_process(rank)+2
          endif
          !Numero de puntos a utilizar en los procesadores


          !Output to the log file           
          if(rank.eq.0) then
             open(200,file=log_file,form='formatted',status='old',position='append')          
             write(200,"()")
             write(200,"(I8A45I8)") grid1_points_process(0)+2 &
                  ,'puntos reales de la malla 1 para el proceso',rank
             
             do i=1,half_numtasks-2
                write(200,"(I8A45I8)") grid1_points_process(i)+2 &
                     ,'puntos reales de la malla 1 para el proceso',i
             enddo
             
             write(200,"(I8A45I8)") grid1_points_process(half_numtasks-1)+1 &
                  ,'puntos reales de la malla 1 para el proceso',half_numtasks-1
             write(200,"()")
             close(200)
          endif
          !Output to the log file           


          !Allocation of memory for first grid
          allocate(x1(1:Nx))
          if(dyn_ec.eq.1) allocate(sgn_x1(1:Nx))
          
          allocate(mod_psi1(1:Nx))
          allocate(mod_pi1(1:Nx))
          
          allocate(auxiliar1(1:Nx))
          allocate(auxiliar2(1:Nx))!se necesita para el metodo 1
          
          allocate(B_mu1(1:Nx))
          allocate(C_R1(1:Nx))
          allocate(C_L1(1:Nx))
          allocate(C_01(1:Nx))
          
          allocate(re_psi1(1:Nx))
          allocate(im_psi1(1:Nx))
          
          allocate(re_pi1(1:Nx))
          allocate(im_pi1(1:Nx))
          
          allocate(re_psi_p1(1:Nx))
          allocate(im_psi_p1(1:Nx))
          allocate(re_pi_p1(1:Nx))
          allocate(im_pi_p1(1:Nx))
          
          allocate(re_psi_a1(1:Nx))
          allocate(im_psi_a1(1:Nx))
          allocate(re_pi_a1(1:Nx))
          allocate(im_pi_a1(1:Nx))
          
          allocate(re_spsi1(1:Nx))
          allocate(im_spsi1(1:Nx))
          allocate(re_spi1(1:Nx))
          allocate(im_spi1(1:Nx))
          
       else if(rank.ge.half_numtasks) then
          
          !Test de Coherencia
          if(epsilon.eq.zero) then
             open(200,file=log_file,form='formatted',status='old',position='append')          
             write(200,"(A70)") 'Error, you do not need two grids when epsilon = zero'
             close(200)
             call MPI_FINALIZE(ierr)
             stop
          endif
          if(epsilon.eq.2.0*mu_0) then
             open(200,file=log_file,form='formatted',status='old',position='append')          
             write(200,"(A70)") 'Error, you do not need two grids when epsilon = 2*mu_0'
             close(200)
             call MPI_FINALIZE(ierr)
             stop
          endif
          !Test de Coherencia


          !Numero de puntos a utilizar en los procesadores
          if(rank.eq.half_numtasks) then
             Nx=grid2_points_process(rank)+2
             
          else if(rank.eq.(numtasks-1)) then
             Nx=grid2_points_process(rank)+1
             
          else
             Nx=grid2_points_process(rank)+2
          endif
          !Numero de puntos a utilizar en los procesadores


          !Output to the log file 
          if(rank.eq.half_numtasks) then
             open(200,file=log_file,form='formatted',status='old',position='append')          
             write(200,"()")
             write(200,"(I8A45I8)") grid2_points_process(half_numtasks)+2 &
                  ,'puntos reales de la malla 1 para el proceso',rank
             
             do i=half_numtasks+1,numtasks-2
                write(200,"(I8A45I8)") grid2_points_process(i)+2 &
                     ,'puntos reales de la malla 1 para el proceso',i
             enddo
             
             write(200,"(I8A45I8)") grid2_points_process(numtasks-1)+1 &
                  ,'puntos reales de la malla 1 para el proceso',numtasks-1
             write(200,"()")
             close(200)
          endif
          !Output to the log file 


          !Allocation of memory for second grid
          allocate(x2(1:Nx))
          if(dyn_ec.eq.1) allocate(sgn_x2(1:Nx))
          
          allocate(mod_psi2(1:Nx))
          allocate(mod_pi2(1:Nx))
          
          allocate(auxiliar1(1:Nx))!se necesita para el metodo 1
          allocate(auxiliar2(1:Nx))
          
          allocate(B_mu2(1:Nx))
          allocate(C_R2(1:Nx))
          allocate(C_L2(1:Nx))
          allocate(C_02(1:Nx))
          
          allocate(re_psi2(1:Nx))
          allocate(im_psi2(1:Nx))
          
          allocate(re_pi2(1:Nx))
          allocate(im_pi2(1:Nx))
          
          allocate(re_psi_p2(1:Nx))
          allocate(im_psi_p2(1:Nx))
          allocate(re_pi_p2(1:Nx))
          allocate(im_pi_p2(1:Nx))
          
          allocate(re_psi_a2(1:Nx))
          allocate(im_psi_a2(1:Nx))
          allocate(re_pi_a2(1:Nx))
          allocate(im_pi_a2(1:Nx))
          
          allocate(re_spsi2(1:Nx))
          allocate(im_spsi2(1:Nx))
          allocate(re_spi2(1:Nx))
          allocate(im_spi2(1:Nx))
          
       endif

    end if

  end subroutine allocate
