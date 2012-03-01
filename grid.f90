! ============================================================================
! grid.f90
! ============================================================================
! Construye la malla para cada proceso

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


  subroutine initial_data

    use arrays
    
    implicit none
    
    include 'mpif.h'

    real(double) :: aux1,aux2

    if(method.eq.1) then

       if(grids_number.eq.1) then
       
          aux1 = one/sqrt(sqrt(cons_pi*two)*delta_mu)
          aux2 = sqrt(three/(16.0*cons_pi*G))*P_phi_star/(mu_star*hbar)
          auxiliar1 = abs(x1)**(three/four)*aux1*exp(-(x1-mu_star)**two /(four*delta_mu**two)) 
          auxiliar2 = abs(x1)**(three/four)*aux1*exp(-(x1+mu_star)**two /(four*delta_mu**two)) 
          
          re_psi1 = auxiliar1*cos(aux2*(x1-mu_star)) &
               + auxiliar2*cos(aux2*(x1+mu_star)) + shift
          
          im_psi1 = -auxiliar1*sin(aux2*(x1-mu_star)) &
               +  auxiliar2*sin(aux2*(x1+mu_star)) + shift
          
          re_pi1 = auxiliar1 * ( sgnp*sqrt(16.0*cons_pi*G/three)*(x1-mu_star)*mu_star & 
               * cos(aux2*(x1-mu_star))/(two*delta_mu**two) + &
               abs(P_phi_star)*sin(aux2*(x1-mu_star))/(hbar*mu_star) ) &
               + auxiliar2 * ( sgnp*sqrt(16.0*cons_pi*G/three)*(-x1-mu_star)*mu_star & 
               * cos(aux2*(-x1-mu_star))/(two*delta_mu**two) + &
               abs(P_phi_star)*sin(aux2*(-x1-mu_star))/(hbar*mu_star) ) + shift
          
          im_pi1 = auxiliar1 * ( abs(P_phi_star)*x1*cos(aux2*(x1-mu_star))/(hbar*mu_star) &
               -sgnp*sqrt(16.0*cons_pi*G/three) & 
               *(x1-mu_star)*mu_star*sin(aux2*(x1-mu_star))/(two*delta_mu**two)) &
               + auxiliar2 * (-abs(P_phi_star)*x1*cos(aux2*(-x1-mu_star))/(hbar*mu_star) &
               -sgnp*sqrt(16.0*cons_pi*G/three) & 
               *(-x1-mu_star)*mu_star*sin(aux2*(-x1-mu_star))/(two*delta_mu**two)) + shift       


       else if(grids_number.eq.2) then

          if(rank.lt.half_numtasks) then

             aux1 = one/sqrt(sqrt(cons_pi*two)*delta_mu)
             aux2 = sqrt(three/(16.0*cons_pi*G))*P_phi_star/(mu_star*hbar)
             auxiliar1 = abs(x1)**(three/four)*aux1*exp(-(x1-mu_star)**two /(four*delta_mu**two)) 
             auxiliar2 = abs(x1)**(three/four)*aux1*exp(-(x1+mu_star)**two /(four*delta_mu**two)) 
             
             re_psi1 = auxiliar1*cos(aux2*(x1-mu_star)) &
                  + auxiliar2*cos(aux2*(x1+mu_star)) + shift
             
             im_psi1 = -auxiliar1*sin(aux2*(x1-mu_star)) &
                  +  auxiliar2*sin(aux2*(x1+mu_star)) + shift
             
             re_pi1 = auxiliar1 * ( sgnp*sqrt(16.0*cons_pi*G/three)*(x1-mu_star)*mu_star & 
                  * cos(aux2*(x1-mu_star))/(two*delta_mu**two) + &
                  abs(P_phi_star)*sin(aux2*(x1-mu_star))/(hbar*mu_star) ) &
                  + auxiliar2 * ( sgnp*sqrt(16.0*cons_pi*G/three)*(-x1-mu_star)*mu_star & 
                  * cos(aux2*(-x1-mu_star))/(two*delta_mu**two) + &
                  abs(P_phi_star)*sin(aux2*(-x1-mu_star))/(hbar*mu_star) ) + shift
             
             im_pi1 = auxiliar1 * ( abs(P_phi_star)*x1*cos(aux2*(x1-mu_star))/(hbar*mu_star) &
                  -sgnp*sqrt(16.0*cons_pi*G/three) & 
                  *(x1-mu_star)*mu_star*sin(aux2*(x1-mu_star))/(two*delta_mu**two)) &
                  + auxiliar2 * (-abs(P_phi_star)*x1*cos(aux2*(-x1-mu_star))/(hbar*mu_star) &
                  -sgnp*sqrt(16.0*cons_pi*G/three) & 
                  *(-x1-mu_star)*mu_star*sin(aux2*(-x1-mu_star))/(two*delta_mu**two)) + shift       
             
          else if(rank.ge.half_numtasks) then
             
             aux1 = one/sqrt(sqrt(cons_pi*two)*delta_mu)
             aux2 = sqrt(three/(16.0*cons_pi*G))*P_phi_star/(mu_star*hbar)
             auxiliar1 = abs(x2)**(three/four)*aux1*exp(-(x2-mu_star)**two /(four*delta_mu**two)) 
             auxiliar2 = abs(x2)**(three/four)*aux1*exp(-(x2+mu_star)**two /(four*delta_mu**two)) 
             
             re_psi2 = auxiliar1*cos(aux2*(x2-mu_star)) &
                  + auxiliar2*cos(aux2*(x2+mu_star)) + shift
             
             im_psi2 = -auxiliar1*sin(aux2*(x2-mu_star)) &
                  +  auxiliar2*sin(aux2*(x2+mu_star)) + shift
             
             re_pi2 = auxiliar1 * ( sgnp*sqrt(16.0*cons_pi*G/three)*(x2-mu_star)*mu_star & 
                  * cos(aux2*(x2-mu_star))/(two*delta_mu**two) +&
                  abs(P_phi_star)*sin(aux2*(x2-mu_star))/(hbar*mu_star) ) &
                  + auxiliar2 * ( sgnp*sqrt(16.0*cons_pi*G/three)*(-x2-mu_star)*mu_star & 
                  * cos(aux2*(-x2-mu_star))/(two*delta_mu**two) +&
                  abs(P_phi_star)*sin(aux2*(-x2-mu_star))/(hbar*mu_star) ) + shift
             
             im_pi2 = auxiliar1 * ( abs(P_phi_star)*x2*cos(aux2*(x2-mu_star))/(hbar*mu_star) &
                  - sgnp*sqrt(16.0*cons_pi*G/three) * & 
                  (x2-mu_star)*mu_star*sin(aux2*(x2-mu_star))/(two*delta_mu**two)) &
                  + auxiliar2 * (-abs(P_phi_star)*x2*cos(aux2*(-x2-mu_star))/(hbar*mu_star) &
                  - sgnp*sqrt(16.0*cons_pi*G/three) * & 
                  (-x2-mu_star)*mu_star*sin(aux2*(-x2-mu_star))/(two*delta_mu**two)) + shift

          end if
          
       endif

    else if(method.eq.2) then       

       if(grids_number.eq.1) then
          
          aux1 = sqrt(three/(16.0*cons_pi*G))*P_phi_star/hbar
          aux2 = sqrt(16.0*cons_pi/three)*((mu_star/delta_mu)**two)/two
          
          auxiliar1 = abs(x1/mu_star)**(one/four) &
               *exp(-0.25*(log(abs(x1/mu_star))*mu_star/delta_mu)**two)
          
          re_psi1 = auxiliar1*cos(-aux1*log(abs(x1/mu_star))) + shift
          im_psi1 = auxiliar1*sin(-aux1*log(abs(x1/mu_star))) + shift
          
          re_pi1  = auxiliar1*(sgnp*aux2*log(abs(x1/mu_star))*cos(-aux1*log(abs(x1/mu_star))) &
               - abs(P_phi_star)*sin(-aux1*log(abs(x1/mu_star)))/hbar ) + shift
          
          im_pi1  = auxiliar1*(sgnp*aux2*log(abs(x1/mu_star))*sin(-aux1*log(abs(x1/mu_star))) &
               + abs(P_phi_star)*cos(-aux1*log(abs(x1/mu_star)))/hbar ) + shift
          
          
          !THIS IS BECAUSE THERE ARE A Log(zero) WHEN epsilon=zero
          if(epsilon.eq.zero) then
             re_psi1(grid_center)=shift!zero
             im_psi1(grid_center)=shift!zero
             
             re_pi1(grid_center)=shift!zero
             im_pi1(grid_center)=shift!zero
          endif
          

       else if(grids_number.eq.2) then

          if(rank.lt.half_numtasks) then
             
             aux1 = sqrt(three/(16.0*cons_pi*G))*P_phi_star/hbar
             aux2 = sqrt(16.0*cons_pi/three)*((mu_star/delta_mu)**two)/two
             
             auxiliar1 = abs(x1/mu_star)**(one/four) &
                  *exp(-0.25*(log(abs(x1/mu_star))*mu_star/delta_mu)**two)
             
             re_psi1 = auxiliar1*cos(-aux1*log(abs(x1/mu_star))) + shift
             im_psi1 = auxiliar1*sin(-aux1*log(abs(x1/mu_star))) + shift
             
             re_pi1  = auxiliar1*(sgnp*aux2*log(abs(x1/mu_star))*cos(-aux1*log(abs(x1/mu_star))) &
                  - abs(P_phi_star)*sin(-aux1*log(abs(x1/mu_star)))/hbar ) + shift

             im_pi1  = auxiliar1*(sgnp*aux2*log(abs(x1/mu_star))*sin(-aux1*log(abs(x1/mu_star))) &
                  + abs(P_phi_star)*cos(-aux1*log(abs(x1/mu_star)))/hbar ) + shift
             
             
             !!Creo que esto no es necesario, no se puede epsilon cero con dos mallas
             !THIS IS BECAUSE THERE ARE A Log(zero) WHEN epsilon=zero
             if(epsilon.eq.zero) then
                re_psi1(grid_center)=shift!zero
                im_psi1(grid_center)=shift!zero
                
                re_pi1(grid_center)=shift!zero
                im_pi1(grid_center)=shift!zero
             endif


          else if(rank.ge.half_numtasks) then
             
             aux1 = sqrt(three/(16.0*cons_pi*G))*P_phi_star/hbar
             aux2 = sqrt(16.0*cons_pi/three)*((mu_star/delta_mu)**two)/two

             auxiliar2 = abs(x2/mu_star)**(one/four) &
                  *exp(-0.25*(log(abs(x2/mu_star))*mu_star/delta_mu)**two)
             
             re_psi2 = auxiliar2*cos(-aux1*log(abs(x2/mu_star))) + shift
             im_psi2 = auxiliar2*sin(-aux1*log(abs(x2/mu_star))) + shift
             
             re_pi2  = auxiliar2*(sgnp*aux2*log(abs(x2/mu_star))*cos(-aux1*log(abs(x2/mu_star))) &
                  - abs(P_phi_star)*sin(-aux1*log(abs(x2/mu_star)))/hbar ) + shift
             im_pi2  = auxiliar2*(sgnp*aux2*log(abs(x2/mu_star))*sin(-aux1*log(abs(x2/mu_star))) &
                  + abs(P_phi_star)*cos(-aux1*log(abs(x2/mu_star)))/hbar ) + shift
             
       endif

    endif

    else if(method.eq.3) then       

       re_psi1 = zero
       re_pi1  = zero
       im_psi1 = zero
       im_pi1  = zero

       if(grids_number.eq.2) then
          re_psi2 = zero
          re_pi2  = zero
          im_psi2 = zero
          im_pi2  = zero
       endif
       
    endif
    
  end subroutine initial_data
  
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


  subroutine grid

    use arrays

    implicit none
    
    include 'mpif.h'

    integer :: i
    real(double) :: aux1,aux2,Xaux1,Xaux2
    integer :: stat(MPI_STATUS_SIZE)
 
    if(dyn_ec.eq.0) then 

       if(grids_number.eq.1) then

          aux1 = cons_pi*G/(9.0*abs(mu_0)**three)
          aux2 = three/two
          
          do i=1,Nx
             x1(i) = grid1_length(rank) + (i-1)*dx*grids_number
          enddo

          !Writing to log file
          if(rank.eq.0) then
             open(200,file=log_file,form='formatted',status='old',position='append')
             
             write(200,"()")
             write(200,"(AI8)") '    punto inicial y final del proceso',rank
             write(200,"(2F15.5)") x1(1),x1(Nx)
             
             do i=1,numtasks-1
                call MPI_RECV(Xaux1,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,stat,ierr) 
                call MPI_RECV(Xaux2,1,MPI_DOUBLE_PRECISION,i,10*i,MPI_COMM_WORLD,stat,ierr) 
                write(200,"()")
                write(200,"(AI8)") '    punto inicial y final del proceso',i
                write(200,"(2F15.5)") Xaux1,Xaux2
             enddo

             close(200)

          else
             call MPI_SEND(x1(1),1,MPI_DOUBLE_PRECISION,0,rank,MPI_COMM_WORLD,ierr)
             call MPI_SEND(x1(Nx),1,MPI_DOUBLE_PRECISION,0,10*rank,MPI_COMM_WORLD,ierr)
          endif
          !Writing to log file
   
          B_mu1 = (2.0/(3.0*mu_0))**6.0 * abs( abs(x1+mu_0)**(3.0/4.0) & 
               -abs(x1-mu_0)**(3.0/4.0) )**6.0
          
          C_R1 = aux1*abs( abs(x1+3.0*mu_0)**aux2 - abs(x1+mu_0)**aux2 ) !C+
          C_L1 = aux1*abs( abs(x1-mu_0)**aux2 - abs(x1-3.0*mu_0)**aux2 )!C-
          C_01 = -C_R1 - C_L1 !C0
       
          
       else if(grids_number.eq.2) then
       
          if(rank.lt.half_numtasks) then

             aux1 = cons_pi*G/(9.0*abs(mu_0)**three)
             aux2 = three/two
              
             do i=1,Nx
                x1(i) = grid1_length(rank) + (i-1)*dx*grids_number 
             enddo

             !Writing to log file
             if(rank.eq.0) then
                open(200,file=log_file,form='formatted',status='old',position='append')

                write(200,"()")
                write(200,"(AI8)") '    punto inicial y final del proceso',rank
                write(200,"(2F15.5)") x1(1),x1(Nx)
                
                do i=1,half_numtasks-1
                   call MPI_RECV(Xaux1,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,stat,ierr) 
                   call MPI_RECV(Xaux2,1,MPI_DOUBLE_PRECISION,i,10*i,MPI_COMM_WORLD,stat,ierr) 
                   write(200,"()")
                   write(200,"(AI8)") '    punto inicial y final del proceso',i
                   write(200,"(2F15.5)") Xaux1,Xaux2
                enddo
                
                close(200)

             else
                call MPI_SEND(x1(1),1,MPI_DOUBLE_PRECISION,0,rank,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x1(Nx),1,MPI_DOUBLE_PRECISION,0,10*rank,MPI_COMM_WORLD,ierr)
             endif
             !Writing to log file

             B_mu1 = (2.0/(3.0*mu_0))**6.0 * abs( abs(x1+mu_0)**(3.0/4.0) & 
                  -abs(x1-mu_0)**(3.0/4.0) )**6.0
             
             C_R1 = aux1*abs( abs(x1+3.0*mu_0)**aux2 - abs(x1+mu_0)**aux2 ) !C+
             C_L1 = aux1*abs( abs(x1-mu_0)**aux2 - abs(x1-3.0*mu_0)**aux2 )!C-
             C_01 = -C_R1 - C_L1 !C0

          
          else  if(rank.ge.half_numtasks) then

             aux1 = cons_pi*G/(9.0*abs(mu_0)**three)
             aux2 = three/two

             do i=1,Nx
                x2(i) = grid2_length(rank) + (i-1)*dx*grids_number 
             enddo
   
             !Writing to log file
             if(rank.eq.half_numtasks) then
                open(200,file=log_file,form='formatted',status='old',position='append')
                
                write(200,"()")
                write(200,"(AI8)") '    punto inicial y final del proceso',rank
                write(200,"(2F15.5)") x2(1),x2(Nx)
                
                do i=half_numtasks+1,numtasks-1
                   call MPI_RECV(Xaux1,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,stat,ierr) 
                   call MPI_RECV(Xaux2,1,MPI_DOUBLE_PRECISION,i,10*i,MPI_COMM_WORLD,stat,ierr) 
                   write(200,"()")
                   write(200,"(AI8)") '    punto inicial y final del proceso',i
                   write(200,"(2F15.5)") Xaux1,Xaux2
                enddo
                
                close(200)

             else
                call MPI_SEND(x2(1),1,MPI_DOUBLE_PRECISION,half_numtasks,rank,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x2(Nx),1,MPI_DOUBLE_PRECISION,half_numtasks,10*rank,MPI_COMM_WORLD,ierr)
             endif
             !Writing to log file
         
             B_mu2 = (2.0/(3.0*mu_0))**6.0 * abs( abs(x2+mu_0)**(3.0/4.0) &
                  -abs(x2-mu_0)**(3.0/4.0) )**6.0
             
             C_R2 = aux1*abs( abs(x2+3.0*mu_0)**aux2 - abs(x2+mu_0)**aux2 ) !C+
             C_L2 = aux1*abs( abs(x2-mu_0)**aux2 - abs(x2-3.0*mu_0)**aux2 )!C-
             C_02 = -C_R2 - C_L2 !C0
          
          endif
       
       end if

    else if (dyn_ec.eq.1) then     

       !NO ESTA LISTO
       
    else if (dyn_ec.eq.2) then     
       
       if(grids_number.eq.1) then
          
          do i=1,Nx
             x1(i) = grid1_length(rank) + (i-1)*dx*grids_number 
          enddo
       
          B_mu1 = x1
          
          C_R1 = x1 + two*lambda !C+
          C_L1 = x1 - two*lambda !C-
          C_01 = -four*x1 !C0
       
       else if(grids_number.eq.2) then
          
          if(rank.lt.half_numtasks) then

             do i=1,Nx
                x1(i) = grid1_length(rank) + (i-1)*dx*grids_number 
             enddo
       
             B_mu1 = x1
             
             C_R1 = x1 + two*lambda !C+
             C_L1 = x1 - two*lambda !C-
             C_01 = -four*x1 !C0
             
          else if(rank.ge.half_numtasks) then

             do i=1,Nx
                x2(i) = grid2_length(rank) + (i-1)*dx*grids_number 
             enddo
             
             B_mu2 = x2
             
             C_R2 = x2 + two*lambda !C+
             C_L2 = x2 - two*lambda !C-
             C_02 = -four*x2 !C0

          endif

       endif

    endif

  end subroutine grid
     
