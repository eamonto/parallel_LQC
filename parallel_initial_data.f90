! ============================================================================
! parallel_initial_data.f90
! ============================================================================
! Aca es donde conocen el rango los procesadores y se dividen los puntos 
! de las mallas entre los procesadores.

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


 subroutine parallel_initial_data

   use arrays
   
   implicit none
   
   include 'mpif.h'
    
   integer aux1,aux2,i

   if(grids_number.eq.1) then
      
      allocate(grid1_points_process(0:numtasks))
      
      aux1=mod(Nx,numtasks) 
      aux2=(Nx-aux1)/numtasks
      
      if(rank.eq.0) then
         open(200,file=log_file,form='formatted',status='old',position='append')
         write(200,"()")
         write(200,"(A20I8)") 'puntos sobrantes ',aux1
         write(200,"(A25I8)") 'puntos por procesador ',aux2
         write(200,"()")
      endif
      
      grid1_points_process=aux2
      
      do i=0,aux1-1
         grid1_points_process(i)=grid1_points_process(i)+1
         if(rank.eq.0)  write(200,"(I8A25I8)") grid1_points_process(i),'puntos para el proceso',i
      end do
      
      do i=numtasks-1,aux1,-1
         if(rank.eq.0)  write(200,"(I8A25I8)") grid1_points_process(i),'puntos para el proceso',i
      end do
      
      if(rank.eq.0) close(200)
      
      
   else if (grids_number.eq.2) then
      
      half_numtasks=numtasks/2
      
      if(rank.eq.0) then
         open(200,file=log_file,form='formatted',status='old',position='append')
         write(200,"(A29I8)"),'mitad de los procesadores',half_numtasks
      endif
      
      !Test de coherencia
      if(mod(numtasks,2).ne.0) then
         if(rank.eq.0) write(200,"(A10)")'--Error--'
         if(rank.eq.0) write(200,"(A40)")'The number of process must be even'
         close(200)
         call MPI_FINALIZE(ierr)
         stop
      end if
      !Test de coherencia
      
      allocate(grid1_points_process(0:half_numtasks))
      allocate(grid2_points_process(half_numtasks:numtasks))
      
      aux1=mod(Nx,half_numtasks)
      aux2=(Nx-aux1)/half_numtasks
      
      if(rank.eq.0) then 
         write(200,"()")
         write(200,"(A21I8)") 'puntos sobrantes ',aux1
         write(200,"(A26I8)") 'puntos por procesador ',aux2
         write(200,"()")
      endif

      grid1_points_process=aux2
      grid2_points_process=aux2
      
      do i=0,aux1-1
         grid1_points_process(i)=grid1_points_process(i)+1
         grid2_points_process(half_numtasks+i)=grid2_points_process(half_numtasks+i)+1
         if(rank.eq.0) then 
            write(200,"(I8A37I8)") grid1_points_process(i),'puntos de la malla 1 para el proceso',i
            write(200,"(I8A37I8)") grid2_points_process(half_numtasks+i) &
                 ,'puntos de la malla 2 para el proceso',half_numtasks+i
         endif
      enddo
      
      if(rank.eq.0) then 
         do i=numtasks-1,half_numtasks+aux1,-1
            write(200,"(I8A37I8)") grid2_points_process(i),'puntos de la malla 2 para el proceso',i
         enddo
         do i=half_numtasks-1,aux1,-1
            write(200,"(I8A37I8)") grid1_points_process(i),'puntos de la malla 1 para el proceso',i
         enddo
         
         write(200,"()")
         close(200)
      endif

   endif

 end subroutine parallel_initial_data
