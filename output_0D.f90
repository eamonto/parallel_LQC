! ============================================================================
! output_0D.f90
! ============================================================================
! Escribe la salida de los observables,
! Solo el proceso cero tiene esta información.

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


  subroutine output_0D

    use arrays

    implicit none

    logical :: firstcall

    data firstcall /.true./
    save firstcall

    if (firstcall) then
       firstcall = .false.
       open(110,file=trim(output_dir)//'/norm.t',form='formatted',status='replace')
       open(111,file=trim(output_dir)//'/volume.t',form='formatted',status='replace')
       open(112,file=trim(output_dir)//'/re_P_phi.t',form='formatted',status='replace')
       open(113,file=trim(output_dir)//'/im_P_phi.t',form='formatted',status='replace')
    else
       open(110,file=trim(output_dir)//'/norm.t',form='formatted',status='old',position='append')
       open(111,file=trim(output_dir)//'/volume.t',form='formatted',status='old',position='append')
       open(112,file=trim(output_dir)//'/re_P_phi.t',form='formatted',status='old',position='append')
       open(113,file=trim(output_dir)//'/im_P_phi.t',form='formatted',status='old',position='append')
    endif

    write(110,"(2ES24.16)") time,norm_sum
    write(111,"(2ES24.16)") time,volume_sum
    write(112,"(2ES24.16)") time,re_P_phi_sum
    write(113,"(2ES24.16)") time,im_P_phi_sum

    close(110)
    close(111)
    close(112)
    close(113)

  end subroutine output_0D
