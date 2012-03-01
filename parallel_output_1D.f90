! ============================================================================
! parallel_output_1D.f90
! ============================================================================
! Coordina la escritura de datos a disco de los procesadores.

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


  subroutine parallel_output_1D

    use arrays

    implicit none

    include 'mpif.h'

    integer :: i,ok,ready
    integer stat(MPI_STATUS_SIZE)

    ok=1
    ready=1

    if(rank.eq.0) then
       call output_1D

       if(numtasks.gt.1) then
          do i=1,numtasks-1
             call MPI_SEND(ready,1,MPI_INTEGER,i,20,MPI_COMM_WORLD,ierr)
             call MPI_RECV(ok,1,MPI_INTEGER,i,26,MPI_COMM_WORLD,stat,ierr) 
         enddo
       end if
    else

       call MPI_RECV(ok,1,MPI_INTEGER,0,20,MPI_COMM_WORLD,stat,ierr)
       call output_1D
       call MPI_SEND(ready,1,MPI_INTEGER,0,26,MPI_COMM_WORLD,ierr)

    endif

  end subroutine parallel_output_1D
