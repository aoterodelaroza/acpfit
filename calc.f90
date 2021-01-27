module calc
  implicit none

  private
  public :: calc_stats
  public :: runfit_inf
  public :: runfit_scanatom
  public :: runfitl
  public :: runfit_manual
  public :: runeval_input
  public :: runeval_file
  public :: runoctavedump
  public :: runoctavedump_binary
  public :: runtest
  private :: choose
  private :: comb

  integer, parameter :: ilong = selected_int_kind(32)

contains
  
  ! calculate the statistics for a given evaluation
  subroutine calc_stats(y,stat,coef,c0)
    use types, only: stats
    use global, only: nfit, w, nset, ytarget, yadd, &
       iset_ini, iset_step, iset_n, naddfiles, addcoef

    real*8, intent(in) :: y(:) ! results of the fit (no dispersion)
    type(stats), intent(inout) :: stat ! statistics of the fit
    real*8, intent(in), optional :: coef(:) ! coefficients
    real*8, intent(in), optional :: c0(:) ! coefficients

    integer :: i, j, id, n
    real*8, allocatable :: yaddl(:), dy(:)

    ! allocate
    allocate(yaddl(nfit),dy(nfit))

    ! calculate the additional energy contribution
    yaddl = 0d0
    do i = 1, naddfiles
       yaddl = yaddl + addcoef(i) * yadd(:,i)
    end do

    ! write down some constants
    if (present(coef)) then
       stat%norm = sqrt(sum((coef*c0)**2))
       stat%norm1 = sum(abs(coef*c0))
       stat%maxcoef = maxval(abs(coef*c0))
    else
       stat%norm = 0d0
       stat%norm1 = 0d0
       stat%maxcoef = 0d0
    end if
    stat%nset = nset

    ! calculate the wrms
    if (abs(sum(w)) < 1d-10) then
       stat%wrms = -1d0
    else
       stat%wrms = sqrt(sum(w * (y + yaddl - ytarget)**2))
    end if

    ! allocate rms and mae
    if (allocated(stat%rms)) deallocate(stat%rms)
    if (allocated(stat%mae)) deallocate(stat%mae)
    if (allocated(stat%mse)) deallocate(stat%mse)
    if (allocated(stat%maewnz)) deallocate(stat%maewnz)
    allocate(stat%rms(nset))
    allocate(stat%mae(nset))
    allocate(stat%mse(nset))
    allocate(stat%maewnz(nset))

    ! calculate rms and mae for all sets
    dy = y + yaddl - ytarget
    do i = 1, nset
       id = iset_ini(i) - iset_step(i)
       stat%rms(i) = 0d0
       stat%mae(i) = 0d0
       stat%mse(i) = 0d0
       stat%maewnz(i) = 0d0

       n = 0
       do j = 1, iset_n(i)
          id = id + iset_step(i)
          stat%rms(i) = stat%rms(i) + (dy(id))**2
          stat%mae(i) = stat%mae(i) + abs(dy(id))
          stat%mse(i) = stat%mse(i) + dy(id)
          if (abs(w(id)) > 1d-14) then
             n = n + 1
             stat%maewnz(i) = stat%maewnz(i) + abs(dy(id))
          end if
       end do
       stat%rms(i) = sqrt(stat%rms(i) / iset_n(i))
       stat%mae(i) = stat%mae(i) / iset_n(i)
       stat%mse(i) = stat%mse(i) / iset_n(i)
       if (n > 0) then
          stat%maewnz(i) = stat%maewnz(i) / n
       else
          stat%maewnz(i) = 0d0
       end if
    end do

  end subroutine calc_stats

  ! Fit the ACP with all possible terms. Good for testing.
  subroutine runfit_inf(outeval,outacp)
    use global, only: ncols, nfitw, nfit, global_printeval, x, global_printacp, col,&
       coef0
    use types, only: stats
    use tools_io, only: uout, faterr, ferror, string
    
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    integer :: idx(ncols), i
    real*8 :: coef(nfitw), c0(nfitw), y(nfit)
    type(stats) :: stat

    ! header
    write (uout,'("+ Generating the ACP that uses all available terms, ACP(Inf)")')

    if (ncols > nfitw) &
       call ferror("runfit_inf","Too many columns ("//string(ncols)//") for current data ("//string(nfitw)//")",faterr)

    ! all possible terms
    idx = 0
    do i = 1, ncols
       idx(i) = i
       c0(i) = coef0(col(i)%iexp,col(i)%l,col(i)%iatom)
    end do

    ! run least squares
    call lsqr(ncols,idx,coef)
    y = matmul(x(:,idx),coef(1:ncols))

    ! print stats and evaluation
    call calc_stats(y,stat,coef(1:ncols),c0(1:ncols))
    call global_printeval("final",y,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",ncols,idx,coef(1:ncols),outacp)

  end subroutine runfit_inf

  ! Fit an ACP using an iterative procedure over the atoms. For each
  ! atom, select the combination of nuse terms that minimizes the
  ! error subject to the maximum 2-norm (maxnorm), maximum 1-norm
  ! (maxnorm1), and/or maximum absolute coefficient (maxcoef)
  ! constraints.
  subroutine runfit_scanatom(nuse,maxnorm,maxnorm1,maxcoef,maxenergy,&
     imaxenergy,outeval,outacp,minl,maxl)
    use global, only: ncols, nfitw, nfit, global_printeval, x, global_printacp, natoms,&
       atom, maxlmax, lmax, nexp, ywtarget, coef0, col, lname
    use types, only: stats, realloc
    use tools_io, only: uout, string, ferror, faterr, string, ioj_left
    
    integer, intent(in) :: nuse
    real*8, intent(in) :: maxnorm
    real*8, intent(in) :: maxnorm1
    real*8, intent(in) :: maxcoef(:,:,:)
    real*8, intent(in) :: maxenergy(:)
    integer, intent(in) :: imaxenergy(:)
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp
    integer, intent(in) :: minl(30), maxl(30)

    character(len=:), allocatable :: str, aux, aux2
    integer :: i, j, k, ncyc, n1, n2, nanuse, nanz, nsame, id
    real*8 :: coef(nfitw), c0(nfitw), y(nfitw), y0(nfit), wrms, minwrms
    real*8 :: norm, norm1, minnorm, minnorm1
    real*8 :: acoef, minacoef, aene(size(imaxenergy,1)), minene(size(imaxenergy,1))
    type(stats) :: stat
    integer, allocatable :: iatidx(:,:), natidx(:), lidx(:,:)
    integer :: co(nuse), idx(natoms*nuse), idx0(natoms*nuse), ll(nuse)
    integer :: idxsave(nuse)
    logical :: cchan, ok

    ! header
    write (uout,'("+ Generating ACP using atom iterations, with ",A," terms per atom")') string(nuse)
    write (uout,*)

    ! some dimensions
    nanuse = natoms*nuse
    if (nanuse > nfitw) &
       call ferror("runfit_inf","Too many columns ("//string(ncols)//") for current data ("//string(nfitw)//")",faterr)

    ! check the maxenergy
    if (imaxenergy(1) /= 0) then
       do i = 1, size(imaxenergy,1)
          if (imaxenergy(i) < 1 .or. imaxenergy(i) > nfit) &
             call ferror("runfit_inf","MAXENERGY system outside trainig set range",faterr)
       end do
    end if

    ! prepare the sets of atomic indices
    allocate(iatidx(maxlmax*nexp,natoms),natidx(natoms),lidx(maxlmax*nexp,natoms))
    natidx = 0
    iatidx = 0
    lidx = 0
    n1 = 0
    do i = 1, natoms
       n2 = 0
       do j = 1, lmax(i)
          do k = 1, nexp
             n1 = n1 + 1
             n2 = n2 + 1
             iatidx(n2,i) = n1
             lidx(n2,i) = j
             natidx(i) = natidx(i) + 1
          end do
       end do
    end do

    ! are we constraining the channels?
    cchan = any(minl /= 0) .or. any(maxl /= huge(1))

    ! run the main loop (over cycles)
    ncyc = 0
    idx0 = 0
    nanz = 0
    nsame = 0
    main: do while(.true.)
       ncyc = ncyc + 1
       write (uout,'("# Cycle ",A)') string(ncyc)

       ! run the loop over atoms
       do i = 1, natoms
          write (uout,'("# Atom ",A," (",A,")")') string(i), string(atom(i))
          idx = idx0
          idxsave = idx0(nuse*(i-1)+1:nuse*i)
          nanz = min(nanz + nuse,nanuse)
          minwrms = huge(1d0)
          minnorm = huge(1d0)
          minnorm1 = huge(1d0)
          minacoef = huge(1d0)
          minene = huge(1d0)
          
          ! Run over all the combinations. To save memory and allow
          ! parallelization, use Buckles' algorithm.
          !$omp parallel do private(coef,c0,y,wrms,norm,norm1,acoef,co,aene,ll,ok) &
          !$omp firstprivate(idx) schedule(dynamic)
          do j = 1, binom(natidx(i),nuse)
             ! get the indices for this combination
             call comb(natidx(i),nuse,j,co)
             if (cchan) then
                ll = lidx(co,i)
                ok = .true.
                do k = 1, lmax(i)
                   ok = ok .and. (count(ll == k) <= maxl(k)) .and. (count(ll == k) >= minl(k))
                end do
                if (.not.ok) cycle
             end if

             ! put them together with the rest of the atoms
             idx(nuse*(i-1)+1:nuse*i) = iatidx(co,i)
             
             ! prepare the c0
             do k = 1, nanz
                c0(k) = coef0(col(idx(k))%iexp,col(idx(k))%l,col(idx(k))%iatom)
             end do

             ! run least squares and calculate wrms
             call lsqr(nanz,idx(1:nanz),coef,y)
             wrms = sqrt(sum((y-ywtarget)**2))
             norm = sqrt(sum((coef(1:nanz) * c0(1:nanz))**2))
             norm1 = sum(abs(coef(1:nanz) * c0(1:nanz)))
             acoef = maxval(abs(coef(1:nanz) * c0(1:nanz)))
             if (imaxenergy(1) > 0) then
                call energy_contrib(nanz,idx(1:nanz),coef(1:nanz),imaxenergy,aene)
             else
                aene = 0d0
             end if

             ! check the maxcoef conditions
             ok = .true.
             do k = 1, nanz
                ok = ok .and. abs(coef(k) * c0(k)) < maxcoef(col(idx(k))%iexp,col(idx(k))%l,col(idx(k))%iatom)
             end do

             ! apply the discard criteria; save minimum wrms
             !$omp critical (save)
             if (ok .and. wrms < minwrms .and. norm < maxnorm .and. norm1 < maxnorm1 .and.&
                all(aene < maxenergy)) then
                idx0 = idx
                minwrms = wrms
                minnorm = norm
                minnorm1 = norm1
                minacoef = acoef
                minene = aene
             end if
             !$omp end critical (save)
          end do
          !$omp end parallel do

          ! check we had at least one combination
          if (minwrms == huge(1d0)) &
             call ferror("runfit_scanatom","Could not find any combination matching the input criteria",faterr)

          ! increase nsame if the new indices are the same as the old ones
          if (all(idxsave == idx0(nuse*(i-1)+1:nuse*i))) then
             nsame = nsame + 1
          else
             nsame = 1
          end if

          ! some output
          do j = 1, natoms
             str = ""
             do k = 1, nuse
                id = idx0(nuse*(j-1)+k)
                if (id > 0) then
                   aux =  string(id,4,ioj_left) // "(" // string(col(id)%atom) // &
                      "," // string(lname(col(id)%l)) // "," // string(col(id)%n) // &
                      "," // trim(string(col(id)%eexp,'f',10,2,ioj_left)) // ")"
                   aux2 = str // "  " // string(aux,20,ioj_left)
                   str = aux2
                end if
             end do
             if (len_trim(str) > 0) &
                write (uout,'(2X,A,2X)') str
          end do
          write (uout,'("  wrms    = ",A)') string(minwrms,'f',14,8)
          write (uout,'("  2-norm  = ",A)') string(minnorm,'f',14,8)
          write (uout,'("  1-norm  = ",A)') string(minnorm1,'f',14,8)
          write (uout,'("  maxcoef = ",A)') string(minacoef,'f',14,8)
          write (uout,'("  maxene  = ",5(A,X))') (string(minene(j),'f',14,8),j=1,size(imaxenergy,1))
          write (uout,'("  nsame   = ",A)') string(nsame)

          ! exit if we're done
          if (nsame == natoms) exit main
       end do
    end do main
    write (uout,*)

    ! final results
    call lsqr(nanuse,idx0,coef)
    y0 = matmul(x(:,idx0),coef(1:nanuse))

    ! prepare the c0
    do k = 1, nanuse
       c0(k) = coef0(col(idx0(k))%iexp,col(idx0(k))%l,col(idx0(k))%iatom)
    end do

    ! print stats and evaluation
    call calc_stats(y0,stat,coef(1:nanuse),c0(1:nanuse))
    call global_printeval("final",y0,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",nanuse,idx0,coef(1:nanuse),outacp)

  end subroutine runfit_scanatom

  ! Fit an ACP using an iterative procedure over the atoms. For each
  ! atom, select the combination of nuse terms that minimizes the
  ! error subject to the maximum 2-norm (maxnorm), maximum 1-norm
  ! (maxnorm1), and/or maximum absolute coefficient (maxcoef)
  ! constraints.
  subroutine runfitl(maxnorm,maxnorm1,maxcoef,maxenergy,imaxenergy,ltop,seq,outeval,outacp)
    use global, only: ncols, nfitw, nfit, global_printeval, x, global_printacp, natoms,&
       atom, maxlmax, lmax, nexp, ywtarget, coef0, col, lname
    use types, only: stats, realloc
    use tools_io, only: uout, string, ferror, faterr, string, ioj_left
    
    real*8, intent(in) :: maxnorm
    real*8, intent(in) :: maxnorm1
    real*8, intent(in) :: maxcoef(:,:,:)
    real*8, intent(in) :: maxenergy(:)
    integer, intent(in) :: imaxenergy(:)
    integer, intent(in) :: ltop(:,:)
    integer, intent(in) :: seq(:,:)
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    integer :: i, j, k, l, m, n1, ncyc, nsame
    real*8 :: coef(nfitw), c0(nfitw), y0(nfit), y(nfitw)
    character(len=:), allocatable :: str, aux, aux2
    real*8 :: wrms, norm, norm1, acoef, aene(size(imaxenergy,1))
    real*8 :: minene(size(imaxenergy,1)), minene2(size(imaxenergy,1))
    real*8 :: minwrms, minnorm, minnorm1, minacoef
    real*8 :: minwrms_, minnorm_, minnorm1_, minacoef_
    type(stats) :: stat
    integer, allocatable :: iatidx(:,:,:), idx0(:,:,:), nidx0(:,:), idx(:)
    integer, allocatable :: co(:), idxroot(:), idxsave(:), idxsave2(:)
    integer :: id, nroot, nn, nord, iord
    logical :: ok, lupdate, saved
    integer, allocatable :: iord_atom(:), iord_chan(:)

    ! header
    write (uout,'("+ Generating ACP using channel + atom iterations")') 
    write (uout,*)

    ! check the maxenergy
    if (imaxenergy(1) /= 0) then
       do i = 1, size(imaxenergy,1)
          if (imaxenergy(i) < 1 .or. imaxenergy(i) > nfit) &
             call ferror("runfit_inf","MAXENERGY system outside training set range",faterr)
       end do
    end if

    ! atom/channel order
    nord = count(seq > 0)
    allocate(iord_atom(nord), iord_chan(nord))
    do i = 1, natoms
       do j = 1, lmax(i)
          iord_atom(seq(j,i)) = i
          iord_chan(seq(j,i)) = j
       end do
    end do
    
    ! prepare the sets of atomic indices
    allocate(iatidx(nexp,maxlmax,natoms))
    iatidx = 0
    n1 = 0
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             n1 = n1 + 1
             iatidx(k,j,i) = n1
          end do
       end do
    end do

    ! prepare the indices for the BSIP
    allocate(idx0(nexp,maxlmax,natoms),nidx0(maxlmax,natoms))
    idx0 = 0
    nidx0 = 0

    ! initialize global quantities
    minwrms = huge(1d0)
    minnorm = huge(1d0)
    minnorm1 = huge(1d0)
    minacoef = huge(1d0)
    minene = huge(1d0)

    ! run the main loop (over cycles)
    ncyc = 0
    nsame = 0
    main: do while(.true.)
       ncyc = ncyc + 1
       write (uout,'("# Cycle ",A)') string(ncyc)

       ! run the loop over atoms
       do iord = 1, nord
          i = iord_atom(iord)
          j = iord_chan(iord)
          write (uout,'("# Atom ",A," (",A,")")') string(i), string(atom(i))
          write (uout,'("# Channel ",A," (",A,")")') string(j), string(lname(j))

          ! build the root index
          nroot = count(idx0 /= 0)
          if (allocated(idxroot)) deallocate(idxroot)
          allocate(idxroot(max(nroot,1)))
          nroot = 0
          do k = 1, natoms
             do l = 1, lmax(k)
                do m = 1, nidx0(l,k)
                   nroot = nroot + 1
                   idxroot(nroot) = idx0(m,l,k)
                end do
             end do
          end do

          ! allocate the save array
          if (allocated(idxsave)) deallocate(idxsave)
          allocate(idxsave(1))
          idxsave = 0

          saved = .false.
          ! run over all possible number of terms for this channel
          do k = 1, min(ltop(j,i),nexp)
             ! allocate the combination array
             if (allocated(co)) deallocate(co)
             allocate(co(k))
             co = 0

             ! allocate the current index array and prepare the slice for the new combinations
             if (allocated(idx)) deallocate(idx)
             allocate(idx(nroot + k))
             nn = nroot + k
             idx(1:nroot) = idxroot
             idx(nroot+1:) = 0

             minwrms_ = huge(1d0)
             minnorm_ = huge(1d0)
             minnorm1_ = huge(1d0)
             minacoef_ = huge(1d0)
             minene2 = huge(1d0)

             ! allocate the second save array
             if (allocated(idxsave2)) deallocate(idxsave2)
             allocate(idxsave2(k))
             idxsave2 = 0

             ! Run over all the combinations with this number of terms
             !$omp parallel do private(coef,c0,y,wrms,norm,norm1,acoef,aene,ok) firstprivate(co,idx)
             do l = 1, binom(nexp,k)
                ! get the indices for this combination
                call comb(nexp,k,l,co)

                ! put them together with the rest of the atoms
                idx(nroot+1:) = 0
                idx(nroot+1:nn) = iatidx(co,j,i)

                ! prepare the c0
                do m = 1, nn
                   c0(m) = coef0(col(idx(m))%iexp,col(idx(m))%l,col(idx(m))%iatom)
                end do

                ! run least squares and calculate wrms
                call lsqr(nn,idx(1:nn),coef,y)
                wrms = sqrt(sum((y-ywtarget)**2))
                norm = sqrt(sum((coef(1:nn) * c0(1:nn))**2))
                norm1 = sum(abs(coef(1:nn) * c0(1:nn)))
                acoef = maxval(abs(coef(1:nn) * c0(1:nn)))
                if (imaxenergy(1) > 0) then
                   call energy_contrib(nn,idx(1:nn),coef(1:nn),imaxenergy,aene)
                else
                   aene = 0d0
                end if

                ! check the maxcoef conditions
                ok = .true.
                do m = 1, nn
                   ok = ok .and. abs(coef(m) * c0(m)) < maxcoef(col(idx(m))%iexp,col(idx(m))%l,col(idx(m))%iatom)
                end do

                ! apply the discard criteria; save minimum wrms
                if (ok .and. wrms < minwrms_ .and. norm < maxnorm .and. &
                   norm1 < maxnorm1 .and. all(aene < maxenergy)) then
                   !$omp critical (save)
                   idxsave2 = idx(nroot+1:nn)
                   minwrms_ = wrms
                   minnorm_ = norm
                   minnorm1_ = norm1
                   minacoef_ = acoef
                   minene2 = aene
                   !$omp end critical (save)
                end if
             end do
             !$omp end parallel do

             ! write down the new combination, if lower than the old one
             if (minwrms_ < minwrms) then
                saved = .true.
                minwrms = minwrms_
                minnorm = minnorm_
                minnorm1 = minnorm1_
                minacoef = minacoef_
                minene = minene2 
                if (allocated(idxsave)) deallocate(idxsave)
                allocate(idxsave(size(idxsave2,1)))
                idxsave = idxsave2
             end if
          end do

          ! process minwrms; only update if we found at least one
          ! successful combination
          if (saved) then
             ! increase nsame if the new indices are the same as the old ones
             k = size(idxsave,1)
             lupdate = (k /= nidx0(j,i))
             if (.not.lupdate) &
                lupdate = any(idx0(1:k,j,i) /= idxsave)

             if (lupdate) then
                nsame = 1
                nidx0(j,i) = k
                idx0(1:k,j,i) = idxsave
             else
                nsame = nsame + 1
             end if
          else
             nsame = nsame + 1
          end if

          ! some output
          do k = 1, natoms
             do l = 1, lmax(k)
                str = ""
                do m = 1, nidx0(l,k)
                   id = idx0(m,l,k)
                   if (id > 0) then
                      aux =  string(id,4,ioj_left) // "(" // string(col(id)%atom) // &
                         "," // string(lname(col(id)%l)) // "," // string(col(id)%n) // &
                         "," // trim(string(col(id)%eexp,'f',10,4,ioj_left)) // ")"
                      aux2 = str // "  " // string(aux,20,ioj_left)
                      str = aux2
                   end if
                end do
                if (len_trim(str) > 0) &
                   write (uout,'(2X,A,2X)') str
             end do
          end do
          write (uout,'("  wrms    = ",A)') string(minwrms,'f',14,8)
          write (uout,'("  2-norm  = ",A)') string(minnorm,'f',14,8)
          write (uout,'("  1-norm  = ",A)') string(minnorm1,'f',14,8)
          write (uout,'("  maxcoef = ",A)') string(minacoef,'f',14,8)
          write (uout,'("  maxene  = ",5(A,X))') (string(minene(k),'f',14,8),k=1,size(imaxenergy,1))
          write (uout,'("  nsame   = ",A)') string(nsame)

          ! exit when we're done
          if (nsame == nord) exit main
          write (uout,*)
       end do
    end do main

    ! build the root index
    if (allocated(idx)) deallocate(idx)
    allocate(idx(count(idx0 /= 0)))
    nn = 0
    do k = 1, natoms
       do l = 1, lmax(k)
          do m = 1, nidx0(l,k)
             nn = nn + 1
             idx(nn) = idx0(m,l,k)
          end do
       end do
    end do

    ! final results
    call lsqr(nn,idx(1:nn),coef)
    y0 = matmul(x(:,idx(1:nn)),coef(1:nn))

    ! prepare the c0
    do m = 1, nn
       c0(m) = coef0(col(idx(m))%iexp,col(idx(m))%l,col(idx(m))%iatom)
    end do

    ! print stats and evaluation
    call calc_stats(y0,stat,coef(1:nn),c0(1:nn))
    call global_printeval("final",y0,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",nn,idx,coef(1:nn),outacp)

  end subroutine runfitl

  ! Find the least-squares ACP with the user-provided ACP terms.
  subroutine runfit_manual(outeval,outacp)
    use global, only: natoms, ncols, col, lname, nfitw, x, nfit, global_printeval,&
       global_printacp, coef0
    use tools_io, only: uin, getline, lgetword, equal, isinteger, isreal, &
       ferror, faterr, lower, string
    use types, only: realloc, stats
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    integer :: n, ninp, lp, i
    integer, allocatable :: idx(:)
    character(len=:), allocatable :: word, line, lstr
    logical :: ok
    real*8 :: rinp, coef(nfitw), c0(nfitw), y(nfit)
    type(stats) :: stat

    ! initialize
    allocate(idx(2*natoms))
    n = 0

    do while (getline(uin,line))
       lp=1
       word = lgetword(line,lp)
       if (equal(word,"end") .or. equal(word,"endrun")) then
          exit
       else
          lstr = lgetword(line,lp)
          ok = isinteger(ninp,line,lp)
          ok = ok .and. isreal(rinp,line,lp)
          if (.not.ok) &
             call ferror("runfit_manual","Invalid ACP term in manual RUN FIT: "//trim(adjustl(line)),faterr)
          
          ok = .false.
          do i = 1, ncols
             if (equal(lower(col(i)%atom),word) .and. &
                equal(lstr,lower(lname(col(i)%l))) .and. col(i)%n == ninp .and.&
                abs(col(i)%eexp - rinp) < 1d-10) then
                n = n + 1
                if (n > size(idx,1)) call realloc(idx,2*n)
                idx(n) = i
                ok = .true.
                exit
             end if
          end do
          if (.not.ok) &
             call ferror("runfit_manual","ACP term not found: "//trim(adjustl(line)),faterr)
       end if
    end do
    call realloc(idx,n)
    
    ! check that we don't have too many terms
    if (n == 0) &
       call ferror("runfit_manual","no ACP terms in manual RUN FIT",faterr)
    if (n > nfitw) &
       call ferror("runfit_manual","Too many columns ("//string(n)//") for current data ("//string(nfitw)//")",faterr)

    ! run least squares
    call lsqr(n,idx,coef)
    y = matmul(x(:,idx),coef(1:n))

    ! prepare the c0
    do i = 1, n
       c0(i) = coef0(col(idx(i))%iexp,col(idx(i))%l,col(idx(i))%iatom)
    end do

    ! print stats and evaluation
    call calc_stats(y,stat,coef(1:n),c0(1:n))
    call global_printeval("final",y,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",n,idx,coef(1:n),outacp)

  end subroutine runfit_manual

  ! Evaluate a given ACP using the current data. From the input file.
  subroutine runeval_input(outeval)
    use global, only: natoms, ncols, col, lname, x, nfit, global_printeval,&
       global_printacp, coef0
    use tools_io, only: uin, getline, lgetword, equal, isinteger, isreal, &
       ferror, faterr, lower, string
    use types, only: realloc, stats
    character*(*), intent(in) :: outeval

    integer :: n, ninp, lp, i
    integer, allocatable :: idx(:)
    character(len=:), allocatable :: word, line, lstr
    logical :: ok
    real*8 :: rinp, y(nfit), r2inp
    real*8, allocatable :: coef(:), c0(:)
    type(stats) :: stat

    ! initialize
    allocate(idx(2*natoms))
    allocate(coef(2*natoms))
    n = 0

    ! parse the input file
    do while (getline(uin,line))
       lp=1
       word = lgetword(line,lp)
       if (equal(word,"end") .or. equal(word,"endrun")) then
          exit
       else
          lstr = lgetword(line,lp)
          ok = isinteger(ninp,line,lp)
          ok = ok .and. isreal(rinp,line,lp)
          ok = ok .and. isreal(r2inp,line,lp)
          if (.not.ok) &
             call ferror("runeval_input","Invalid term+coef in RUN EVAL: "//trim(adjustl(line)),faterr)
          
          ok = .false.
          do i = 1, ncols
             if (equal(lower(col(i)%atom),word) .and. &
                equal(lstr,lower(lname(col(i)%l))) .and. col(i)%n == ninp .and.&
                abs(col(i)%eexp - rinp) < 1d-10) then
                n = n + 1
                if (n > size(idx,1)) then
                   call realloc(idx,2*n)
                   call realloc(coef,2*n)
                end if
                idx(n) = i
                coef(n) = r2inp / coef0(col(i)%iexp,col(i)%l,col(i)%iatom)
                ok = .true.
                exit
             end if
          end do
          if (.not.ok) &
             call ferror("runeval_input","ACP term not found: "//trim(adjustl(line)),faterr)
       end if
    end do
    call realloc(idx,n)
    call realloc(coef,n)
    
    ! check that we don't have too many terms
    if (n == 0) &
       call ferror("runeval_input","no ACP terms in manual RUN EVAL",faterr)

    ! evaluate
    y = matmul(x(:,idx),coef(1:n))

    ! prepare the c0
    allocate(c0(n))
    do i = 1, n
       c0(i) = coef0(col(idx(i))%iexp,col(idx(i))%l,col(idx(i))%iatom)
    end do

    ! print stats and evaluation
    call calc_stats(y,stat,coef(1:n),c0(1:n))
    call global_printeval("final",y,stat,outeval)
    
  end subroutine runeval_input

  ! Evaluate a given ACP using the current data. From an external file.
  subroutine runeval_file(outeval,inacp)
    use global, only: natoms, lname, x, nfit, global_printeval,&
       global_printacp, coef0, whichatom, whichcol, whichexp, col, naddfiles,&
       addfile, addcoef
    use tools_io, only: getline, lgetword, equal, isinteger, isreal, &
       ferror, faterr, lower, string, fopen_read, fclose
    use types, only: realloc, stats
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: inacp

    integer :: nblock
    integer :: n, lp, j, lu, iatom, iang, nterms
    integer, allocatable :: idx(:)
    character(len=:), allocatable :: word, line
    logical :: ok, found
    real*8 :: y(nfit)
    real*8, allocatable :: coef(:), c0(:)
    type(stats) :: stat
    integer :: idum, i
    real*8 :: rdum1, rdum2

    ! initialize
    allocate(idx(2*natoms))
    allocate(coef(2*natoms))
    n = 0

    ! parse the external file
    lu = fopen_read(inacp)
    do while (getline(lu,line))
       lp=1
       
       ! Read the additional contribution coefficients. Else, skip comments
       if (index(line,"Additional energy contributions") > 0) then
          do while (getline(lu,line))
             lp = index(line,"!")
             if (lp > 0) then
                lp = lp + 1
                word = lgetword(line,lp)
                lp = index(word,":")
                word = trim(word(1:lp-1))

                found = .false.
                do i = 1, naddfiles
                   if (equal(word,trim(addfile(i)))) then
                      lp = index(line,"=") + 1
                      ok = isreal(addcoef(i),line,lp)
                      found = .true.
                   end if
                end do
                if (.not.found) exit
             end if
          end do
          lp = 1
       end if

       ! skip the comments
       if (line(lp:lp) == "!") cycle

       ! read the atom
       word = lgetword(line,lp)
       iatom = whichatom(word)
       if (iatom == 0) &
          call ferror("runeval_file","Unknown atom: "//trim(word),faterr)
       
       ! read the number of blocks
       ok = getline(lu,line,.true.)
       lp = 1
       word = lgetword(line,lp)
       ok  = isinteger(nblock,line,lp)
       if (.not.ok) &
          call ferror("runeval_file","Could not parse line: "//trim(line),faterr)
       nblock = nblock + 1

       ! read the blocks
       do iang = 1, nblock
          ! label
          ok = getline(lu,line,.true.)
          lp = 1
          word = lgetword(line,lp)
          if (lname(iang) /= word(1:1)) &
             call ferror("runeval_file","Unknown ang. mom. label: "//trim(line),faterr)

          ! number of terms
          ok = getline(lu,line,.true.)
          lp = 1
          ok = isinteger(nterms,line,lp)
          if (.not.ok) &
             call ferror("runeval_file","Incorrect number of terms: "//trim(line),faterr)

          ! read the terms
          do j = 1, nterms
             n = n + 1
             if (n > size(idx,1)) then
                call realloc(idx,2*n)
                call realloc(coef,2*n)
             end if
             ok = getline(lu,line,.true.)
             lp = 1
             ok = isinteger(idum,line,lp)
             ok = ok .and. isreal(rdum1,line,lp)
             ok = ok .and. isreal(rdum2,line,lp)
             if (.not.ok) &
                call ferror("runeval_file","Incorrect term: "//trim(line),faterr)

             idx(n) = whichcol(iatom,iang,idum,rdum1)
             if (idx(n) == 0) &
                call ferror("runeval_file","No data found for term: "//trim(line),faterr)
             coef(n) = rdum2 / coef0(whichexp(idum,rdum1),iang,iatom)
          end do
       end do
    end do
    call fclose(lu)
    call realloc(idx,n)
    call realloc(coef,n)
    
    ! check that we don't have too many terms
    if (n == 0) &
       call ferror("runeval_file","no ACP terms in manual RUN EVAL acp.file",faterr)

    ! evaluate
    y = matmul(x(:,idx),coef(1:n))

    ! prepare the c0
    allocate(c0(n))
    do i = 1, n
       c0(i) = coef0(col(idx(i))%iexp,col(idx(i))%l,col(idx(i))%iatom)
    end do

    ! print stats and evaluation
    call calc_stats(y,stat,coef(1:n),c0(1:n))
    call global_printeval("final",y,stat,outeval)
    
  end subroutine runeval_file

  ! Write the octave dump file for the LASSO fitting script
  subroutine runoctavedump(imode,maxcoef)
    use global, only: natoms, atom, lmax, nexp, eexp, xw, ywtarget, ywadd, coef0,&
       imode_octavedump, imode_octavedump_universal, imode_octavedump_universal_local,&
       addfile
    use tools_io, only: uout, faterr, ferror, string, fopen_write, fclose
    integer, intent(in) :: imode
    real*8, intent(in) :: maxcoef(:,:,:)

    character(len=:), allocatable :: ofile 
    integer :: lu, lmaxx
    integer :: i, j, k, n, m
    real*8, allocatable :: coef0map(:), xwaux(:,:)

    ofile = "octavedump.m"
    write (uout,'("+ Writing octave dump file to octavedump.m"/)')

    lu = fopen_write(ofile)

    if (imode == imode_octavedump) then
       write (lu,'("atoms={...")')
       do i = 1, natoms
          write (lu,'("""",A,""",...")') trim(atom(i))
       end do
       write (lu,'("};")')
    else
       write (lu,'("atoms={""X""};")')
    end if

    if (imode == imode_octavedump) then
       write (lu,'("lmax=[...")')
       write (lu,'(4X,999(A,X))') (string(lmax(i)),i=1,natoms)
       write (lu,'("];")')
    elseif (imode == imode_octavedump_universal) then
       lmaxx = minval(lmax(1:natoms))
       write (lu,'("lmax=[",A,"];")') string(lmaxx)
    elseif (imode == imode_octavedump_universal_local) then
       lmaxx = 1
       write (lu,'("lmax=[",A,"];")') string(lmaxx)
    end if
    write (lu,'("lname={""l"",""s"",""p"",""d"",""f"",""g""};")')

    write (lu,'("explist=[...")')
    write (lu,'(4X,999(A,X))') (trim(string(eexp(i),'f',20,10)),i=1,nexp)
    write (lu,'("];")')

    write (lu,'("nrows = ",A,";")') string(size(xw,1))

    if (imode == imode_octavedump) then
       write (lu,'("ncols = ",A,";")') string(size(xw,2))
    else
       write (lu,'("ncols = ",A,";")') string(lmaxx * nexp)
    end if

    if (any(maxcoef /= huge(1d0))) then
       if (imode == imode_octavedump_universal .or. imode == imode_octavedump_universal_local) &
          call ferror("runoctavedump","maxcoef not compatible with universal/universal_local",faterr)
       write (lu,'("maxcoef=[...")')
       do i = 1, natoms
          do j = 1, lmax(i)
             do k = 1, nexp
                write (lu,'(1X,A)') trim(string(maxcoef(k,j,i),'e',25,14))
             end do
          end do
       end do
       write (lu,'("]'';")')
    end if

    allocate(coef0map(size(xw,2)))
    n = 0
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             n = n + 1
             coef0map(n) = coef0(k,j,i)
          end do
       end do
    end do

    write (lu,'("x=[...")')
    if (imode == imode_octavedump) then
       do i = 1, size(xw,1)
          write (lu,'(1X,99999(A,X))') (trim(string(xw(i,j)/coef0map(j),'e',25,14)),j=1,size(xw,2))
       end do
    else
       allocate(xwaux(size(xw,1),lmaxx * nexp))
       xwaux = 0d0
       n = 0
       do i = 1, natoms
          m = 0
          do j = 1, lmax(i)
             do k = 1, nexp
                n = n + 1
                m = m + 1
                if (j > lmaxx) cycle
                xwaux(:,m) = xwaux(:,m) + xw(:,n) / coef0map(n)
             end do
          end do
       end do
       
       do i = 1, size(xwaux,1)
          write (lu,'(1X,99999(A,X))') (trim(string(xwaux(i,j),'e',25,14)),j=1,size(xwaux,2))
       end do
       deallocate(xwaux)
    end if
    write (lu,'("];")')
    deallocate(coef0map)

    write (lu,'("y=[...")')
    do i = 1, size(ywtarget)
       write (lu,'(1X,A)') trim(string(ywtarget(i),'e',25,14))
    end do
    write (lu,'("];")')

    if (allocated(ywadd)) then
       write (lu,'("yaddnames={...")')
       do j = 1, size(ywadd,2)
          write (lu,'(2X,"""",A,""",...")') trim(addfile(j))
       end do
       write (lu,'("};")')

       write (lu,'("yadd=[...")')
       do i = 1, size(ywadd,1)
          write (lu,'(1X,999(A,X))') (trim(string(ywadd(i,j),'e',25,14)),j=1,size(ywadd,2))
       end do
       write (lu,'("];")')
    end if

    call fclose(lu)

  end subroutine runoctavedump

  ! Write the octave dump file for the LASSO fitting script (binary version)
  subroutine runoctavedump_binary(maxcoef)
    use global, only: natoms, atom, lmax, nexp, eexp, x, ywtarget, yadd, coef0,&
       addfile, w, yref, yempty, ydisp
    use tools_io, only: uout, faterr, ferror, string, fopen_write, fclose
    real*8, intent(in) :: maxcoef(:,:,:)

    character(len=:), allocatable :: ofile , methodname
    integer :: lu, lmaxx
    integer :: i, j, k, n, m
    real*8, allocatable :: coef0map(:), xwaux(:,:), x_(:,:)
    integer*8 :: nrows, ncols, nadd, addmaxl, sizes(7), nmaxc
    character*2, allocatable :: atnames(:)

    ofile = "octavedump.dat"
    write (uout,'("+ Writing binary octave dump file to octavedump.dat"/)')
    lu = fopen_write(ofile,form='unformatted',stream=.true.)

    ! sizes
    nadd = 0
    if (allocated(yadd)) nadd = size(yadd,2)
    nrows = size(x,1)
    ncols = size(x,2)
    addmaxl = 0
    do i = 1, size(yadd,2)
       addmaxl = max(addmaxl,len_trim(addfile(i)))
    end do
    sizes = (/int(natoms,8),int(nexp,8),nrows,ncols,nadd+1,nadd,addmaxl/)
    write (lu) sizes

    ! atomic names
    write (lu) atom(1:natoms)

    ! additional method names
    allocate(character(len=addmaxl)::methodname)
    do i = 1, nadd
       methodname = trim(addfile(i))
       methodname = adjustl(methodname)
       write (lu) methodname
    end do
    
    ! lmax
    write (lu) int(lmax,1)

    ! exponents
    write (lu) eexp(1:nexp)

    ! weights
    write (lu) w

    ! acp terms
    allocate(coef0map(size(x,2)))
    n = 0
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             n = n + 1
             coef0map(n) = coef0(k,j,i)
          end do
       end do
    end do
    allocate(x_(size(x,1),size(x,2)))
    do i = 1, size(x,2)
       x_(:,i) = x(:,i) / coef0map(i)
    end do
    write (lu) x_
    deallocate(x_,coef0map)

    ! reference energy
    write (lu) yref

    ! empty energy
    write (lu) yempty

    ! added energies
    if (nadd > 0) then
       write (lu) yadd
    end if

    ! subtracted energies
    write (lu) ydisp

    if (any(maxcoef /= huge(1d0))) then
       write (lu) int(1,1)
       nmaxc = 0
       do i = 1, natoms
          do j = 1, lmax(i)
             do k = 1, nexp
                nmaxc = nmaxc + 1
             end do
          end do
       end do
       write (lu) nmaxc
       do i = 1, natoms
          do j = 1, lmax(i)
             do k = 1, nexp
                write (lu) maxcoef(k,j,i)
             end do
          end do
       end do
    else
       write (lu) int(0,1)
    end if

    call fclose(lu)

  end subroutine runoctavedump_binary

  ! Write a tetsting ACP that uses all terms and has coefficients that roughtly
  ! give the same average contribution to the wrms over the whole set.
  subroutine runtest(outeval,outacp)
    use global, only: ncols, nfitw, global_printeval, global_printacp, w, yempty, &
       yref, x, nfit, coef0, col
    use types, only: stats
    use tools_io, only: uout, faterr, ferror, string
    
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    integer :: idx(ncols), i
    real*8 :: coef(ncols), c0(ncols), wmae0, y(nfit)
    type(stats) :: stat

    ! header
    write (uout,'("+ Generating testing ACP with all available terms")')

    if (ncols > nfitw) &
       call ferror("runtest","Too many columns ("//string(ncols)//") for current data ("//string(nfitw)//")",faterr)

    ! all possible terms
    idx = 0
    do i = 1, ncols
       idx(i) = i
    end do

    ! empty/ref wmae
    wmae0 = sum(w * abs(yempty-yref)) / sum(w)
    write (uout,'("  wmae(empty/ref) = ",A/)') string(wmae0,'f',12,6)
    wmae0 = wmae0 / ncols * 10d0

    ! find the coefficients
    do i = 1, ncols
       coef(i) = wmae0 / (sum(w * abs(x(:,i))) / sum(w))
       c0(i) = coef0(col(i)%iexp,col(i)%l,col(i)%iatom)
    end do

    ! run least squares
    y = matmul(x(:,idx),coef(1:ncols))

    ! print stats and evaluation
    call calc_stats(y,stat,coef(1:ncols),c0(1:ncols))
    call global_printeval("final",y,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",ncols,idx,coef(1:ncols),outacp)

  end subroutine runtest

  ! Run least squares using ndim columns given by idx. Returns the
  ! coefficients in coef(1:ndim).
  subroutine lsqr(ndim,idx,coef,yout)
    use global, only: natoms, ncols, nfitw, xw, ywtarget
    integer, intent(in) :: ndim
    integer, intent(in) :: idx(ndim)
    real*8, intent(out) :: coef(nfitw)
    real*8, intent(out), optional :: yout(nfitw)

    ! work data for lapack
    integer :: jpvt(natoms*ncols) 
    real*8 :: work(natoms*ncols + 3 * nfitw + 1), xwork(nfitw,ndim)
    integer :: lwork, rank, info

    ! initialize lapack
    jpvt = 0
    lwork = natoms*ncols + 3 * nfitw + 1

    ! build the work slice
    ! xwork(:,1:ndim) = xw(:,idx)
    xwork = xw(:,idx)
    coef = ywtarget

    ! run the least squares 
    call dgelsy(nfitw,ndim,1,xwork,nfitw,coef,nfitw,jpvt,1d-10,rank,work,lwork,info)
    if (present(yout)) then
       xwork = xw(:,idx)
       yout = matmul(xwork(:,1:ndim),coef(1:ndim))
    end if

  end subroutine lsqr

  !> Calculate the energy terms for the ACP given by indices idx with
  !> coefficients coef (each of those with dimension ndim) and return
  !> the maximum energy contribution from all the terms in absolute
  !> value.
  subroutine energy_contrib(ndim,idx,coef,imaxenergy,aene)
    use global, only: nfit, x
    integer, intent(in) :: ndim
    integer, intent(in) :: idx(ndim)
    real*8, intent(in) :: coef(ndim)
    integer, intent(in) :: imaxenergy(:)
    real*8, intent(out) :: aene(size(imaxenergy,1))

    integer :: i
    
    aene = 0d0
    do i = 1, ndim
       aene = max(aene,abs(x(imaxenergy,idx(i))*coef(i)))
    end do

  end subroutine energy_contrib

  !> Calculate n choose k using the multiplicative formula.
  function choose(n,k)
    integer*8 :: choose
    integer, intent(in) :: n, k
 
    integer, parameter :: qp = selected_real_kind(32)
    integer :: i, s
    real(kind=qp) :: x
 
    if (k > n .or. k < 0) then
       choose = 0
       return
    end if

    s = min(k,n-k)
    x = 1.
    do i = 1, s
       x = (x * (n + 1 - i)) / i
    end do
    choose = nint(x,8)

  end function choose
  
  !> Given a combination of k in n elements (idx), calculate the next
  !> combination in the lexicographical order.
  subroutine nextcomb(idx,n,k)
    integer, intent(in) :: n, k
    integer, intent(inout) :: idx(k)

    integer :: i, j
    
    do i = k, 1, -1
       if (idx(i) < n-k+i)  then
          idx(i) = idx(i) + 1
          do j = i+1, k
             idx(j) = idx(j-1)+1
          end do
          return
       end if
    end do

  end subroutine nextcomb

  ! From TOMS/515
  ! Generates a vector from a lexicographical index That is, let C1,
  ! C2, ... Cm be the set of combinations of n items taken p at a
  ! time arranged in lexographical order. Given an integer i, this
  ! routine finds Ci.
  ! B.P. Buckles and M. Lybanon, ACM TOMS 3 (1977) 180-182
  pure subroutine comb(N, P, L, C)                                       
    ! THIS SUBROUTINE FINDS THE COMBINATION SET OF N THINGS
    ! TAKEN P AT A TIME FOR A GIVEN LEXICOGRAPHICAL INDEX.
    ! N - NUMBER OF THINGS IN THE SET
    ! P - NUMBER OF THINGS IN EACH COMBINATION
    ! L - LEXICOGRAPHICAL INDEX OF COMBINATION SOUGHT
    ! C - OUTPUT ARRAY CONTAINING THE COMBINATION SET
    ! THE FOLLOWING RELATIONSHIPS MUST EXIST AMONG THE INPUT
    ! VARIABLES.  L MUST BE GREATER THAN OR EQUAL TO 1 AND LESS
    ! THAN OR EQUAL TO THE MAXIMUM LEXICOGRAPHICAL INDEX.
    ! P MUST BE LESS THAN OR EQUAL TO N AND GREATER THAN ZERO.
    integer, intent(in) :: n, p, l
    integer, intent(out) :: c(p)

    integer :: K, R, P1, I

    ! SPECIAL CASE CODE IF P = 1
    IF(P.EQ.1)THEN
       C(1) = L
       RETURN
    ENDIF

    ! INITIALIZE LOWER BOUND INDEX AT ZERO
    K = 0
    ! LOOP TO SELECT ELEMENTS IN ASCENDING ORDER
    P1 = P - 1
    C(1) = 0
    DO I=1,P1
       ! SET LOWER BOUND ELEMENT NUMBER FOR NEXT ELEMENT VALUE
       IF (I.NE.1) C(I) = C(I-1)
       ! LOOP TO CHECK VALIDITY OF EACH ELEMENT VALUE
10     C(I) = C(I) + 1
       R = BINOM(N-C(I),P-I)
       K = K + R
       IF (K.LT.L) GO TO 10
       K = K - R
    end DO
    C(P) = C(P1) + L - K

  end subroutine comb

  pure function binom(M, N)
    ! ACM ALGORITHM 160 TRANSLATED TO FORTRAN.  CALCULATES THE
    ! NUMBER OF COMBINATIONS OF M THINGS TAKEN N AT A TIME.
    integer, intent(in) :: m, n
    integer :: binom

    INTEGER :: P, I, N1, R

    N1 = N
    P = M - N1
    IF (N1.GE.P) GO TO 10
    P = N1
    N1 = M - P
10  R = N1 + 1
    IF (P.EQ.0) R = 1
    IF (P.LT.2) GO TO 30
    DO I=2,P
       R = (R*(N1+I))/I
    end DO
30  BINOM = R
  end FUNCTION BINOM

end module calc
