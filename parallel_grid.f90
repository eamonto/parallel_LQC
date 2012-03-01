! ============================================================================
! parallel_grid.f90
! ============================================================================
! Halla las longitudes de cada parte de la malla para cada procesador

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


 subroutine parallel_grid

   use arrays
   
   implicit none
   
   include 'mpif.h'

   integer :: i
   
   if(rank.eq.0) then
      open(200,file=log_file,form='formatted',status='old',position='append')
      write(200,"(A32F15.5)") 'longitud total de la malla  ',length
      close(200)
   endif
   
   if(grids_number.eq.1) then
      
      allocate(grid1_length(0:numtasks))
       
      do i=0,numtasks-1
         grid1_length(i)=-length/two +epsilon + i*length/numtasks
      end do
      
   else if(grids_number.eq.2) then
      
      allocate(grid1_length(0:half_numtasks))
      allocate(grid2_length(half_numtasks:numtasks))
      
      do i=0,half_numtasks-1
         !epsilon is necessary because de grid is no centered
         grid1_length(i)=-length/two +epsilon + i*length/half_numtasks
         grid2_length(half_numtasks+i)=-length/two -epsilon + i*length/half_numtasks
      end do
      
   end if
   
 end subroutine parallel_grid
