! ============================================================================
! main.f90
! ============================================================================
! Programa principal

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


  program main

    use param
    use arrays

    implicit none
   
    include 'mpif.h'

    logical :: internal
    integer :: l,k,Ntime
    real(double) :: flag_observables

    flag_observables=zero
 
    !REVISAR ESTA CONDICION 
    !Necessary condition for make the grid symmetric
    if((dyn_ec.eq.0).and.(epsilon.eq.2.0*mu_0)) Nx=Nx-1
    if((dyn_ec.eq.1).and.(epsilon.eq.2.0*mu_0)) Nx=Nx-1
    if((dyn_ec.eq.2).and.(epsilon.eq.2.0*lambda)) Nx=Nx-1
    !REVISAR ESTA CONDICION 

    if (dyn_ec.eq.0) then
       dx = four*mu_0/grids_number     
    else if (dyn_ec.eq.1) then
       dx = four*mu_0/grids_number 
    else if (dyn_ec.eq.2) then
       dx = four*lambda/grids_number
    else
       print *,'--Error-- The parameter of dynamical equations is wrong'
       stop
    endif
   
    length = dx*grids_number*Nx

    !MPI Initialization
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
    !MPI Initialization

    !Log file creation
    log_file='log_'//trim(output_dir)
    if(rank.eq.0) then
       open(200,file=log_file,form='formatted',status='replace')
       close(200)
    endif
    !Log file creation

    if(rank.eq.0) then
       open(200,file=log_file,form='formatted',status='old',position='append')
       write(200,"()")
       write(200,"(A)") '    =========================================='
       write(200,"(A)") '    Starting To Solve the Big Bang Singularity'
       write(200,"(A)") '    =========================================='
       write(200,"()")
       close(200)

       if(stdout.eq.1) then
          write(*,"()")
          write(*,"(A)") ' =========================================='
          write(*,"(A)") ' Starting To Solve the Big Bang Singularity'
          write(*,"(A)") ' =========================================='
          write(*,"()")
       endif

    endif


    !!Aca es donde conocen el rango los procesadores
    !!Se dividen los puntos de las mallas entre los procesadores
    call parallel_initial_data

    if(rank.eq.0) then
       call system('mkdir -p '//trim(output_dir))
       call system('cp param.f90 '//trim(output_dir))
    endif

    !!Alocacion de la memoria para cada procesador
    call allocate

    !!Inicializa todos las variables a cero o valores predefinidos
    call initialize

    !!Halla las longitudes de cada parte de la malla
    !!para cada procesador
    call parallel_grid

    !****Nota: las subrutinas parallel_grid y grid se puede unir

    !!Construye la malla para cada proceso
    call grid

    !!Le da los valores iniciales a las funciones de onda
    !!en cada parte de la malla (o cada proceso)
    call initial_data

    !!Calcula el valor de los oservables: volume,P_phi (y la norma)
    !!cada proceso calcula su parte y luego se unifican los valores
    call observables

    !!Escribe la salida de los observables,
    !!Solo el proceso cero tiene esta informacion
    if(rank.eq.0) call output_0D

    !!Escribe el output de la funcion de onda y la derivada
    call parallel_output_1D


    if(grids_number.eq.1) then
       if(numtasks.gt.1) internal=.true.

       if(rank.eq.0) then 
          open(200,file=log_file,form='formatted',status='old',position='append')
          write(200,"()")
          write(200,"(A)")'    Internal communication ON'    
          close(200)
       endif

    else if(grids_number.eq.2) then
       if(numtasks.ge.4) internal=.true.

       if(rank.eq.0) then 
          open(200,file=log_file,form='formatted',status='old',position='append')
          write(200,"()")
          write(200,"(A)")'    Internal communication ON'    
          close(200)
       endif

    else
       internal=.false.

       if(rank.eq.0) then 
          open(200,file=log_file,form='formatted',status='old',position='append')
          write(200,"()")
          write(200,"(A)")'    Internal communication OFF'    
          close(200)
       endif
    endif

    if(rank.eq.0) then
       open(200,file=log_file,form='formatted',status='old',position='append')
       write(200,"()")
       write(200,"(A)")' -------------------------------'
       write(200,"(A)")' |  Time step  |     Time      |'
       write(200,"(A)")' -------------------------------'
       write(200,"(A5,I6,A6,ES12.5,A3)") ' |   ',0,'    | ',time,'  |'
       close(200)

       if(stdout.eq.1) then
          write(*,"()")
          write(*,"(A)")' -------------------------------'
          write(*,"(A)")' |  Time step  |     Time      |'
          write(*,"(A)")' -------------------------------'
          write(*,"(A5,I6,A6,ES12.5,A3)") ' |   ',0,'    | ',time,'  |'
       endif

    endif


    Ntime=int(Total_time/dt)+1 

     do l=1,Ntime

       time = time + direction*dt
       
       if(rank.eq.0) then
          if(mod(l,every_1D).eq.0) then
             open(200,file=log_file,form='formatted',status='old',position='append')
             write(200,"(A5,I6,A6,ES12.5,A3)") ' |   ',l,'    | ',time,'  |'
             close(200)
             if(stdout.eq.1) write(*,"(A5,I6,A6,ES12.5,A3)") ' |   ',l,'    | ',time,'  |'
          endif
       endif

       call store_levels

       do k=1,4
          
          call sources
 
          if (internal) call internal_boundaries !comunicacion entre las fronteras de los nodos
         
          call external_boundaries !frontera externa total
         
          call evolution(k)

       enddo
       
       if(mod(l,every_0D).eq.0) then
          call observables
          if(rank.eq.0) call output_0D
          flag_observables=one
       endif

       if(mod(l,every_1D).eq.0) then
          if(flag_observables.eq.zero) call observables
          call parallel_output_1D
       endif

       flag_observables=zero

    enddo

    if(rank.eq.0) then
       open(200,file=log_file,form='formatted',status='old',position='append')
       write(200,"(A)")' ------------------------------- '
       write(200,"()")
       write(200,"(A)")' The Big Bang Singularity is Solved '
       write(200,"()")
       close(200)

       if(stdout.eq.1) then
          write(*,"(A)")' ------------------------------- '
          write(*,"()")
          write(*,"(A)")' The Big Bang Singularity is Solved, Congratulations! '
          write(*,"()")
       endif

    endif

    call MPI_FINALIZE(ierr)
    
  end program main
