! This program performs anatomy calculation based on Sparse Structure Learning (Ide et. al. 2008)
! for Givin Two Data sets of Time Serie

program a_ssl
  use standardlize
  use graphical_Lasso

  implicit none

  integer :: i, j, k, l, m
  integer :: ndim   ! (m)  # of dimension
  integer :: nframea, nframeb ! (n)  # of frame

  double precision,allocatable, dimension(:,:) :: xa, xb    ! (m x n), time series
  double precision,allocatable, dimension(:)   :: avxa,avxb ! (m), average of data

  double precision,allocatable,dimension(:,:) :: Aa, Ab, Unit    ! (m x m), presicion matrix
  double precision,allocatable,dimension(:,:) :: iAa, iAb  ! (m x m), inverse presis matrix
  double precision,allocatable,dimension(:,:) :: Sa, Sb    ! (m x m), sample variance-covariance matrix
  double precision,allocatable,dimension(:)   :: an_ab, an_ba        ! (m), anomaly for each dimension
  double precision r ! degree of sparsity
  double precision, allocatable, dimension(:,:) :: B(:,:)   ! dummy matrix
  double precision v(5) ! dummy values
  
  character dataa_filename*200, datab_filename*200
  ! character out_filename*50

  character(10) :: argv, programname
  integer :: argc, iargc

  NAMELIST/assl/ndim, r !, aflag
  open(16,file='parameters_assl')
  read(16,assl)
  close(16)

  argc=iargc()
  call getarg(0,programname)
  if (argc < 2) then !  if (argc < 3) then
     call usage(programname)
  else
     call getarg(1,dataa_filename)
     call getarg(2,datab_filename)
     !     call getarg(3,out_filename)
  end if

  open(21,file=dataa_filename,status='old')
  allocate(xa(ndim,1))
  call read_dat(21, nframea, ndim, xa)

  open(21,file=datab_filename,status='old')
  allocate(xb(ndim,1))
  call read_dat(21, nframeb, ndim, xb)

  ! open(19,file=out_filename,status='new')
  write(*,'("       Number of the data set A = ", I4)'), nframea
  write(*,'("       Number of the data set B = ", I4)'), nframeb
  write(*,'("       Degrees of freedom = ", I4)'), ndim
  write(*,'("       Value of sparsity parameter is ",F8.3)'), r
  
  allocate(avxa(ndim), avxb(ndim))
  allocate(Sa(ndim,ndim), Sb(ndim,ndim))

  allocate(Aa(ndim,ndim), Ab(ndim,ndim))
  allocate(iAa(ndim,ndim), iAb(ndim,ndim))

  call standradlize_ts(xa,nframea,ndim)      ! standradlize the data
  call ave_timeseries(xa,avxa,nframea,ndim)   ! avx <- x
  call vcv_timeseries(xa,nframea,ndim,avxa,Sa) ! S <- x

  call standradlize_ts(xb,nframeb,ndim)      ! standradlize the data
  call ave_timeseries(xb,avxb,nframeb,ndim)   ! avx <- x
  call vcv_timeseries(xb,nframeb,ndim,avxb,Sb) ! S <- x

  deallocate(xa, xb)
  deallocate(avxa,avxb)
  
  k = 0
  write(*,'("The sample variance-covarucance matrix of the data set A")')
  do i=1,ndim,1
     do j=i,ndim,1 !do j=1,ndim,1
        k = k + 1
        write(*,'(F8.3,1X)', advance='no'), Sa(i,j)
        if (mod(k,5) .eq. 0) then
           write(*,'(/)',advance='no')
        end if
     end do
  end do
  write(*,'(/)',advance='no')

  k = 0
  write(*,'("The sample variance-covarucance matrix of the data set B")')
  do i=1,ndim,1
     do j=i,ndim,1 !do j=1,ndim,1
        k = k + 1
        write(*,'(F8.3,1X)', advance='no'), Sb(i,j)
        if (mod(k,5) .eq. 0) then
           write(*,'(/)',advance='no')
        end if
     end do
  end do
  write(*,'(/)',advance='no')
  
  call block_coordinate_descent_method(Aa,iAa,Sa,r,ndim) ! A, iA, <- S, rou, m

  call block_coordinate_descent_method(Ab,iAb,Sb,r,ndim) ! A, iA, <- S, rou, m
  
  k = 0
  write (*,'("The MAP estimation of presicion matrix from the data set A")')
  do i=1,ndim
     do j=i,ndim !do j=1,ndim
        k = k + 1
        write (*,'(F8.3,1X)',advance='no'), Aa(i,j)
        if (mod(k,5) .eq. 0) then
           write(*,'(/)',advance='no')
        end if
     end do
  end do
  write(*,'(/)',advance='no')
  k = 0
  write (*,'("The inverse matrix of the presicion matrix from the data set A")')
  do i=1,ndim
     do j=i,ndim !do j=1,ndim
        k = k + 1
        write (*,'(F8.3,1X)',advance='no'), iAa(i,j)
        if (mod(k,5) .eq. 0) then
           write(*,'(/)',advance='no')
        end if
     end do
  end do
  write(*,'(/)',advance='no')

  k = 0
  write (*,'("The MAP estimation of presicion matrix from the data set B")')
  do i=1,ndim
     do j=i,ndim !do j=1,ndim
        k = k + 1
        write (*,'(F8.3,1X)',advance='no'), Ab(i,j)
        if (mod(k,5) .eq. 0) then
           write(*,'(/)',advance='no')
        end if
     end do
  end do
  write(*,'(/)',advance='no')
  k = 0
  write (*,'("The inverse matrix of the presicion matrix from the data set B")')
  do i=1,ndim
     do j=i,ndim !do j=1,ndim
        k = k + 1
        write (*,'(F8.3,1X)',advance='no'), iAb(i,j)
        if (mod(k,5) .eq. 0) then
           write(*,'(/)',advance='no')
        end if
     end do
  end do
  write(*,'(/)',advance='no')

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocate(Unit(ndim,ndim))                           ! ! !
  !                                                     ! ! !
  ! Unit = 0.0d0                                        ! ! !
  ! do i = 1,ndim                                       ! ! !
  !    do j = 1,ndim                                    ! ! !
  !       do k = 1,ndim                                 ! ! !
  !          Unit(i,j) = Unit(i,j) + Aa(i,k)*iAa(k,j)   ! ! !
  !       end do                                        ! ! !
  !    end do                                           ! ! !
  ! end do                                              ! ! !
  !                                                     ! ! !
  ! do i=1,ndim                                         ! ! !
  !    do j=1,ndim                                      ! ! !
  !       write (*,'(F8.3,1X)',advance='no'), Unit(i,j) ! ! !
  !    end do                                           ! ! !
  !    write(*,'(/)',advance='no')                      ! ! !
  ! end do                                              ! ! !
  ! write(*,'(/)',advance='no')                         ! ! !
  !                                                     ! ! !
  ! deallocate(Unit)                                    ! ! !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  allocate(an_ab(ndim), an_ba(ndim))

!  call calc_anomaly(Aa,Ab,Sa,an_ab,ndim) ! an_ab <- Aa, Ab, Sa, m
  call calc_anomaly2(Aa,Ab,iAa,iAb,Sa,an_ab,ndim) ! an_ab <- Aa, Ab, Sa, m

!  call calc_anomaly(Ab,Aa,Sb,an_ba,ndim) ! an_ba <- Ab, Aa, Sb, m
  call calc_anomaly2(Ab,Aa,iAb,iAa,Sb,an_ba,ndim) ! an_ba <- Ab, Aa, Sb, m
  
  deallocate(Sa, Sb)
  deallocate(Aa, Ab)
  deallocate(iAa, iAb)

  write (*,'("The anomaly of each dimension")')
  write (*,'("  #      a->b     b->a     max)")')
  do i=1,ndim,1
     if (an_ab(i) > an_ba(i)) then
        write(*,'(I3,1X,F8.3,1X,F8.3,1X,F8.3)'),  i,an_ab(i),an_ba(i),an_ab(i)
     else
        write(*,'(I3,1X,F8.3,1X,F8.3,1X,F8.3)'),  i,an_ab(i),an_ba(i),an_ba(i)
     end if
  end do
  ! close(19)

  deallocate(an_ab,an_ba)
  
contains

  ! Show help mesage
  subroutine usage(programname)
    implicit none
    character(10) :: programname
    !    print *,"Usage: ",programname," dataa_filename datab_filename out_filename"
    print *,"Usage: ",programname," dataa_filename datab_filename"
    stop
  end subroutine usage
  
  ! Read input file
  subroutine read_dat (unitnum, N, ndim, x)

    implicit none
    integer, intent(in) :: unitnum
    integer, intent(out) :: N
    integer, intent(in) :: ndim
    double precision,allocatable, dimension(:,:), intent(inout) :: x

    integer :: i, j, k
    double precision, allocatable, dimension(:,:) :: B
    double precision v
    integer err
    
    i = 1
    j = 0
    do 
       read (unitnum, '(F8.3)', advance="NO",iostat=err, end=999), v
       if (err .eq. 0 ) then
          j = j + 1
          if ( j > ndim) then
             allocate(B(ndim,i))
             B = x
             i = i + 1
             j = 1
             deallocate(x)
             allocate(x(ndim,i))
             x(1:ndim,1:i-1) = B(1:ndim,1:i-1)
             deallocate(B)
             x(j,i) = v
          else
             x(j,i) = v
          end if
       end if
    end do
999 N = i - 1

!  do i = 1,Nndx,1
!     write(*,'(I4,1X)', advance="NO"), atom_list(i)
!     if (mod(i,10) .eq. 0) then
!        write(*,"()")
!     end if
!  end do
    
  end subroutine read_dat

end program a_ssl
