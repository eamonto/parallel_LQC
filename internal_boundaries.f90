! ============================================================================
! internal_boundaries.f90
! ============================================================================
! Comunica los valores de las fronteras la malla a los vecinos de cada 
! proceso.

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


subroutine internal_boundaries

  use arrays
  
  implicit none
  
  include 'mpif.h'
  
  integer :: stat(MPI_STATUS_SIZE)

     
  if(grids_number.eq.1) then

     !!Envio de datos hacia la derecha
     if(rank.eq.0) then
        call MPI_SEND(re_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(im_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+2,MPI_COMM_WORLD,ierr)
           
     else if (rank.lt.(numtasks-1)) then
        call MPI_RECV(re_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+1,MPI_COMM_WORLD,stat,ierr) 
        call MPI_RECV(im_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+2,MPI_COMM_WORLD,stat,ierr) 

        call MPI_SEND(re_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(im_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+2,MPI_COMM_WORLD,ierr)

     else if (rank.eq.(numtasks-1)) then
        call MPI_RECV(re_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+1,MPI_COMM_WORLD,stat,ierr) 
        call MPI_RECV(im_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+2,MPI_COMM_WORLD,stat,ierr) 
           
     endif
     
     !!Envio de datos hacia la izquierda
     if(rank.eq.(numtasks-1)) then
        call MPI_SEND(re_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(im_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+2,MPI_COMM_WORLD,ierr)
        
     else if (rank.gt.0) then
        call MPI_RECV(re_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+1,MPI_COMM_WORLD,stat,ierr) 
        call MPI_RECV(im_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+2,MPI_COMM_WORLD,stat,ierr) 

        call MPI_SEND(re_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(im_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+2,MPI_COMM_WORLD,ierr)

     else if (rank.eq.0) then
        call MPI_RECV(re_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+1,MPI_COMM_WORLD,stat,ierr) 
        call MPI_RECV(im_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+2,MPI_COMM_WORLD,stat,ierr) 
        
     endif
        

  else if(grids_number.eq.2) then 
     
     !!Envio de datos hacia la derecha
     if(rank.lt.half_numtasks) then
        
        if(rank.eq.0) then
           call MPI_SEND(re_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+2,MPI_COMM_WORLD,ierr)
              
        else if (rank.lt.(half_numtasks-1)) then
           call MPI_RECV(re_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+1,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+2,MPI_COMM_WORLD,stat,ierr) 
           
           call MPI_SEND(re_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi1(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+2,MPI_COMM_WORLD,ierr)
              
        else if (rank.eq.(half_numtasks-1)) then
           call MPI_RECV(re_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+1,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi1(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+2,MPI_COMM_WORLD,stat,ierr) 
              
        endif
     
     !!Envio de datos hacia la derecha
     else  if(rank.ge.half_numtasks) then
           
        if(rank.eq.half_numtasks) then
           call MPI_SEND(re_spi2(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi2(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+2,MPI_COMM_WORLD,ierr)
           
        else if (rank.lt.(numtasks-1)) then
           call MPI_RECV(re_spi2(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+1,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi2(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+2,MPI_COMM_WORLD,stat,ierr) 
           
           call MPI_SEND(re_spi2(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi2(Nx-1),1,MPI_DOUBLE_PRECISION,rank+1,100+10*rank+2,MPI_COMM_WORLD,ierr)
           
        else if (rank.eq.(numtasks-1)) then
           call MPI_RECV(re_spi2(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+1,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi2(1),1,MPI_DOUBLE_PRECISION,rank-1,100+10*(rank-1)+2,MPI_COMM_WORLD,stat,ierr) 
           
        endif
     endif


     !!Envio de datos hacia la izquierda
     if(rank.lt.half_numtasks) then
           
        if(rank.eq.(half_numtasks-1)) then
           call MPI_SEND(re_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+3,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+4,MPI_COMM_WORLD,ierr)
           
        else if (rank.gt.0) then
           call MPI_RECV(re_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+3,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+4,MPI_COMM_WORLD,stat,ierr) 

           call MPI_SEND(re_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+3,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi1(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+4,MPI_COMM_WORLD,ierr)

        else if (rank.eq.0) then
           call MPI_RECV(re_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+3,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi1(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+4,MPI_COMM_WORLD,stat,ierr) 
           
        endif
        
     !!Envio de datos hacia la izquierda           
     else  if(rank.ge.half_numtasks) then
           
        if(rank.eq.(numtasks-1)) then
           call MPI_SEND(re_spi2(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+3,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi2(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+4,MPI_COMM_WORLD,ierr)
           
        else if (rank.gt.half_numtasks) then
           call MPI_RECV(re_spi2(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+3,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi2(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+4,MPI_COMM_WORLD,stat,ierr) 
           
           call MPI_SEND(re_spi2(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+3,MPI_COMM_WORLD,ierr)
           call MPI_SEND(im_spi2(2),1,MPI_DOUBLE_PRECISION,rank-1,100+10*rank+4,MPI_COMM_WORLD,ierr)

        else if (rank.eq.half_numtasks) then
           call MPI_RECV(re_spi2(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+3,MPI_COMM_WORLD,stat,ierr) 
           call MPI_RECV(im_spi2(Nx),1,MPI_DOUBLE_PRECISION,rank+1,100+10*(rank+1)+4,MPI_COMM_WORLD,stat,ierr) 
           
        endif
     endif
  endif

end subroutine internal_boundaries
