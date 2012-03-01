! ============================================================================
! output_1D.f90
! ============================================================================
! Escribe el output de la función de onda y la derivada de la función de onda.

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


  subroutine output_1D

    use arrays

    implicit none

    include 'mpif.h'

    integer :: i
    logical :: firstcall

    data firstcall /.true./
    save firstcall

    if(rank.eq.0) then
    
       if (firstcall) then
          firstcall = .false.
          open(101,file=trim(output_dir)//'/mod_psi.x',form='formatted',status='replace')
          open(102,file=trim(output_dir)//'/mod_pi.x',form='formatted',status='replace')
       else
          open(101,file=trim(output_dir)//'/mod_psi.x',form='formatted',status='old',position='append')
          open(102,file=trim(output_dir)//'/mod_pi.x',form='formatted',status='old',position='append')
       endif
       
    else
       open(101,file=trim(output_dir)//'/mod_psi.x',form='formatted',status='old',position='append')
       open(102,file=trim(output_dir)//'/mod_pi.x',form='formatted',status='old',position='append')
    endif

    if(rank.eq.0) then
       write(101,"(A8,ES24.16)") '#Time = ',time
       write(102,"(A8,ES24.16)") '#Time = ',time
    endif
    
    !aca se estan escribiendo dos veces el punto 1 y Nx de cada malla
    !A excepcion de los puntos externos de toda la malla
    !No creo que sea gran problema
    do i=1,Nx
       if(mod((i-1),every_1D_grid).eq.0) then
          
          if(grids_number.eq.1) then

             if(mod_psi1(i).lt.shift) mod_psi1(i) = shift
             if(mod_pi1(i).lt.shift)  mod_pi1(i)  = shift
             write(101,"(2ES24.16)") x1(i),mod_psi1(i)
             write(102,"(2ES24.16)") x1(i),mod_pi1(i)

          else if(grids_number.eq.2) then
             
             if(rank.lt.half_numtasks) then

                if(mod_psi1(i).lt.shift) mod_psi1(i) = shift
                if(mod_pi1(i).lt.shift)  mod_pi1(i)  = shift
                write(101,"(2ES24.16)") x1(i),mod_psi1(i)
                write(102,"(2ES24.16)") x1(i),mod_pi1(i)
                
             else if(rank.ge.half_numtasks) then
                
                if(mod_psi2(i).lt.shift) mod_psi2(i) = shift
                if(mod_pi2(i).lt.shift)  mod_pi2(i)  = shift
                write(101,"(2ES24.16)") x2(i),mod_psi2(i)
                write(102,"(2ES24.16)") x2(i),mod_pi2(i)
             
             endif

          endif

       endif

    enddo
    
    if(rank.eq.(numtasks-1)) then
       write(101,*)
       write(102,*)
       
       write(101,*)
       write(102,*)
    end if

    close(101)
    close(102)

  end subroutine output_1D
