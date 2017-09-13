! Copyright (c) 2017 Alberto Otero de la Roza <aoterodelaroza@gmail.com>
!
! acpfit is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! acpfit is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.
program acpfit
  use calc, only: calc_stats, runfit_inf, runfit_scanatom, runfitl, runfit_manual, &
     runeval_input, runeval_file, runtest, runoctavedump
  use files, only: makefilenames, readfiles
  use global, only: fileroot, nset, iset_label, iset_ini, iset_n, iset_step, nfit,&
     outempty, outeval, outacp, imode_fit, imode_fitl, imode_fit_manual, imode_eval,&
     imode_eval_file, imode_test, imode_octavedump, global_init, global_input, global_fillcol,&
     global_printinfo, global_printeval, inacp
  use types, only: realloc, stats
  use tools_io, only: equal, lgetword, getline, uout, isreal, isinteger, getword, &
     string, nwarns, ncomms, ioinit, stdargs, help_me, tictac, &
     start_clock, print_clock
  implicit none

  character(len=:), allocatable :: optv
  type(stats) :: statempty
  real*8, allocatable :: ydum(:), fit_maxenergy(:), fit_maxcoef(:,:,:)
  integer :: ifit_n, imode, maxl(30), minl(30)
  real*8 :: fit_maxnorm, fit_maxnorm1
  integer, allocatable :: imaxenergy(:), ltop(:,:), seq(:,:)

  ! energy terms files in manual input
  integer :: nefilesi ! number of files
  character*255, allocatable :: efilei(:) ! name of files

  ! open input and output files
  call ioinit()
  call stdargs(optv,fileroot)

  ! help message and banner
  write (uout,'("* ACPFIT: fit atom-centered potentials using least squares")')
  if (index(optv,'h') > 0) then
     call help_me()
     goto 999
  end if
  
  ! timing information
  call tictac('ACPFIT started')
  call start_clock()
  write (uout,*)

  ! initialize variables
  call global_init()

  ! read and parse the input
  call global_input(nefilesi,efilei,imode,ifit_n,fit_maxnorm,fit_maxnorm1,fit_maxcoef,&
     fit_maxenergy,imaxenergy,minl,maxl,ltop,seq)

  ! build the file names and process user input re specific energy term file names
  call makefilenames(nefilesi,efilei)

  ! fill the column information and use the noexp to prune the exponent list
  call global_fillcol(fit_maxcoef)

  ! add the "all" set
  nset = nset + 1
  call realloc(iset_label,nset)
  call realloc(iset_ini,nset)
  call realloc(iset_n,nset)
  call realloc(iset_step,nset)
  iset_label(nset) = "all"
  iset_ini(nset) = 1
  iset_n(nset) = nfit
  iset_step(nset) = 1

  ! write the input data information
  call global_printinfo(fit_maxnorm,fit_maxnorm1,fit_maxcoef,fit_maxenergy,imaxenergy)

  ! read the information from the external files
  call readfiles()

  ! calculate the statistics for the empty
  allocate(ydum(nfit))
  ydum = 0d0
  call calc_stats(ydum,statempty)
  call global_printeval("empty",ydum,statempty,outempty)

  ! run the calculation
  if (imode == imode_fit) then
     if (ifit_n < 0) then
        ! all terms in the ACP
        call runfit_inf(outeval,outacp)
     elseif (ifit_n > 0) then
        call runfit_scanatom(ifit_n,fit_maxnorm,fit_maxnorm1,fit_maxcoef,fit_maxenergy,&
           imaxenergy,outeval,outacp,minl,maxl)
     end if
  elseif (imode == imode_fitl) then
     call runfitl(fit_maxnorm,fit_maxnorm1,fit_maxcoef,fit_maxenergy,imaxenergy,&
        ltop,seq,outeval,outacp)
  elseif (imode == imode_fit_manual) then
     call runfit_manual(outeval,outacp)
  elseif (imode == imode_eval) then
     call runeval_input(outeval)
  elseif (imode == imode_eval_file) then
     call runeval_file(outeval,inacp)
  elseif (imode == imode_test) then
     call runtest(outeval,outacp)
  elseif (imode == imode_octavedump) then
     call runoctavedump()
  end if

  write (uout,'("ACPFIT ended succesfully (",A," WARNINGS, ",A," COMMENTS)")')&
     string(nwarns), string(ncomms)
  call print_clock()
  call tictac('ACPFIT finished')

999 continue

end program acpfit

