module calc
  implicit none

  private
  public :: calc_stats
  public :: runfit_inf

contains
  
  ! calculate the statistics for a given evaluation
  subroutine calc_stats(y,stat,coef)
    use types, only: stats
    use global, only: nfit, maxnamelen, w, ydisp, yref, names, nset, ytarget, &
       iset_ini, iset_step, iset_n, coef0

    real*8, intent(in) :: y(:) ! results of the fit (no dispersion)
    type(stats), intent(inout) :: stat ! statistics of the fit
    real*8, intent(in), optional :: coef(:) ! coefficients

    integer :: i, j, id
    real*8 :: dy(nfit)

    ! write down some constants
    if (present(coef)) then
       stat%norm = sqrt(sum(coef**2)) * coef0
       stat%maxcoef = maxval(abs(coef)) * coef0
    else
       stat%norm = 0d0
       stat%maxcoef = 0d0
    end if
    stat%nset = nset

    ! calculate the wrms
    if (abs(sum(w)) < 1d-10) then
       stat%wrms = -1d0
    else
       stat%wrms = sqrt(sum(w * (y-ytarget)**2) / sum(w))
    end if

    ! allocate rms and mae
    if (allocated(stat%rms)) deallocate(stat%rms)
    if (allocated(stat%mae)) deallocate(stat%mae)
    allocate(stat%rms(nset))
    allocate(stat%mae(nset))

    ! calculate rms and mae for all sets
    dy = y - ytarget
    do i = 1, nset
       id = iset_ini(i) - iset_step(i)
       stat%rms(i) = 0d0
       stat%mae(i) = 0d0
       do j = 1, iset_n(i)
          id = id + iset_step(i)
          stat%rms(i) = stat%rms(i) + (dy(id))**2
          stat%mae(i) = stat%mae(i) + abs(dy(id))
       end do
       stat%rms(i) = sqrt(stat%rms(i) / iset_n(i))
       stat%mae(i) = stat%mae(i) / iset_n(i)
    end do

  end subroutine calc_stats

  ! Fit the ACP with all possible terms. Good for testing.
  subroutine runfit_inf(outeval,outacp)
    use global, only: ncols, nfitw, nfit, global_printeval, x, yempty, global_printacp
    use types, only: stats
    
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    integer :: idx(ncols), i
    real*8 :: coef(nfitw), y(nfit)
    type(stats) :: stat

    ! all possible terms
    idx = 0
    do i = 1, ncols
       idx(i) = i
    end do

    ! run least squares
    call lsqr(ncols,idx,coef)
    y = matmul(x(:,idx),coef)

    ! print stats and evaluation
    call calc_stats(y,stat,coef)
    call global_printeval("final",y,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",ncols,idx,coef,outacp)

  end subroutine runfit_inf

  ! Run least squares using ndim columns given by idx. Returns the
  ! coefficients in coef(1:ndim).
  subroutine lsqr(ndim,idx,coef)
    use global, only: natoms, ncols, nfitw, xwork, xw, ywtarget
    integer, intent(in) :: ndim
    integer, intent(in) :: idx(ndim)
    real*8, intent(out) :: coef(nfitw)

    ! work data for lapack
    integer :: jpvt(natoms*ncols) 
    real*8 :: work(natoms*ncols + 3 * nfitw + 1) 
    integer :: lwork, rank, info

    ! initialize lapack
    jpvt = 0
    lwork = natoms*ncols + 3 * nfitw + 1

    ! build the work slice
    xwork(:,1:ndim) = xw(:,idx)
    coef = ywtarget

    ! run the least squares 
    call dgelsy(nfitw,ndim,1,xwork(:,1:ndim),nfitw,coef,nfitw,jpvt,1d-10,rank,work,lwork,info)

  end subroutine lsqr

end module calc
