! This program performs anatomy calculation based on Sparse Structure Learning (Ide et. al. 2008)
! for Givin Two Data sets of Time Serie

program a_ssl
  use standardlize
  use graphical_Lasso

  implicit none

  integer :: i, j, k, l, m

  integer :: ndim   ! (m)  # of dimension
  integer :: nframea, nframeb ! (n)  # of frame

  integer :: ninia, ninib, nf
  integer :: nfina, nfinb
  
  double precision,allocatable, dimension(:,:) :: xa, xb    ! (m x n), time series
  double precision,allocatable, dimension(:)   :: avxa,avxb ! (m), average of data

  double precision,allocatable,dimension(:,:) :: Aa, Ab, Unit   ! (m x m), presicion matrix
  double precision,allocatable,dimension(:,:) :: iAa, iAb       ! (m x m), inverse presis matrix
  double precision,allocatable,dimension(:,:) :: Sa, Sb         ! (m x m), sample variance-covariance matrix
  double precision,allocatable,dimension(:)   :: an_ab, an_ba   ! (m), anomaly for each dimension
  double precision r ! degree of sparsity
  double precision, allocatable, dimension(:,:) :: B            ! dummy matrix
  double precision v(5) ! dummy values
  
  character dataa_filename*200, datab_filename*200
  ! character out_filename*50

  integer :: mode   ! if mdoe = 0, standarization will be done, if mode = 1, not.
  integer :: dbmode ! if mdoe = 0, default, if mode = 1, debug mode.
  
  character(30) :: argv, programname
  character*500 :: arg(20), buff
  integer :: argc, iargc

  !  NAMELIST/assl/ndim, r !, aflag
  !  open(16,file='parameters_assl')
  !  read(16,assl)
  !  close(16)

  ! Default value
  mode = 0
  dbmode = 0
  ninia = 1
  ninib = 1
  nfina = -1
  nfinb = -1
  ndim = 164
  r = 0.30
  ! Options
  argc=iargc()
  call getarg(0,programname)
  if (argc < 2) then !  if (argc < 3) then
     call usage(programname)
  else
     do i=1,argc-2
        call getarg(i,arg(i))
        buff=arg(i)
        if (buff(1:5) .eq. "-mode") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),mode;
           cycle;
        end if
        if (buff(1:7) .eq. "-dbmode") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),mode;
           cycle;
        end if
        if (buff(1:6) .eq. "-ninia") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),ninia;
           cycle;
        end if
        if (buff(1:6) .eq. "-ninib") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),ninib;
           cycle;
        end if
        if (buff(1:6) .eq. "-nfina") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),nfina;
           cycle;
        end if
        if (buff(1:6) .eq. "-nfinb") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),nfinb;
           cycle;
        end if
        if (buff(1:9) .eq. "-Sparsity") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),r;
           cycle;
        end if
        if (buff(1:5) .eq. "-ndim") then
           call getarg(i+1,arg(i+1))
           read(arg(i+1),*),ndim;
           cycle;
        end if
     end do

     call getarg(argc-1, dataa_filename)
     call getarg(argc  , datab_filename)
     !     call getarg(3,out_filename)
  end if
  
  print *,"       Program Name:",programname
  print *,"       Data A is ",dataa_filename
  print *,"       Data B is ",datab_filename

  print *,"       "
  
  open(21,file=dataa_filename,status='old')
  allocate(xa(ndim,1))
  call read_dat(21, nframea, ndim, xa)
  close(21)

  write(*,'  ("       Number of the data set A = ", I4)'), nframea

  if ( nfina .eq. -1) then
     nfina = nframea
  end if

!  print *, nframea, nfina, ninia

  if ( ninia <= 0 ) then
     print *, "Initial step must be larger than 0."
     stop
  end if
  
  if ( nfina < ninia) then
     print *, "Final step must be larger than initial step."
     stop
  end if
  
  nf = nfina - ninia + 1
  allocate(B(ndim,nf))
  B(1:ndim,1:nf) = xa(1:ndim,ninia:nfina)
  deallocate(xa)
  nframea = nf
  allocate(xa(ndim,nframea))
  xa = B
  deallocate(B)

  open(22,file=datab_filename,status='old')
  allocate(xb(ndim,1))
  call read_dat(22, nframeb, ndim, xb)
  close(22)

  write(*,  '("       Number of the data set B = ", I4)'), nframeb

  if ( ninib <= 0 ) then
     print *, "Initial step must be larger than 0."
     stop
  end if
  
  if ( nfinb .eq. -1) then
     nfinb = nframeb
  end if

!  print *, nframeb, nfinb, ninib

  if ( nfinb < ninib) then
     print *, "Final step must be larger than initial step."
     stop
  end if
  
  nf = nfinb - ninib + 1
  allocate(B(ndim,nf))
  B(1:ndim,1:nf) = xb(1:ndim,ninib:nfinb)
  deallocate(xb)
  nframeb = nf
  allocate(xb(ndim,nframeb))
  xb = B
  deallocate(B)

  print  *,  "       Options "
  print  *,  "       Mode                        = ", mode
  print  *,  "       Initial step of data A      = ", ninia
  print  *,  "       Initial step of data B      = ", ninib
  print  *,  "       Final   step of data A      = ", nfina
  print  *,  "       Final   step of data B      = ", nfinb
  print  *,  "       Degrees of freedom          = ", ndim
  write( *,'("        Value of sparsity parameter = ",F8.3)'), r

  print *,"       "
    
  allocate(avxa(ndim), avxb(ndim))
  allocate(Sa(ndim,ndim), Sb(ndim,ndim))

  allocate(Aa(ndim,ndim), Ab(ndim,ndim))
  allocate(iAa(ndim,ndim), iAb(ndim,ndim))

  call standradlize_ts(xa,nframea,ndim,mode)   ! standradlize the data
  call ave_timeseries(xa,avxa,nframea,ndim)    ! avx <- x
  call vcv_timeseries(xa,nframea,ndim,avxa,Sa) ! S <- x

  call standradlize_ts(xb,nframeb,ndim,mode)   ! standradlize the data
  call ave_timeseries(xb,avxb,nframeb,ndim)    ! avx <- x
  call vcv_timeseries(xb,nframeb,ndim,avxb,Sb) ! S <- x

  deallocate(xa, xb)
  deallocate(avxa,avxb)
  
  if ( dbmode .eq. 1 ) then                                        
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
        do j=i,ndim,1                                                       
           k = k + 1                                                        
           write(*,'(F8.3,1X)', advance='no'), Sb(i,j)                      
           if (mod(k,5) .eq. 0) then                                        
              write(*,'(/)',advance='no')                                   
           end if
        end do
     end do
     write(*,'(/)',advance='no')
  end if
  
  write(*,'("The block coordinate descent method will be started.")')

  call block_coordinate_descent_method(Aa,iAa,Sa,r,ndim) ! A, iA, <- S, rou, m

  write(*,'("The block coordinate descent method will be started.")')

  call block_coordinate_descent_method(Ab,iAb,Sb,r,ndim) ! A, iA, <- S, rou, m
  
  if (dbmode .eq. 1) then                                        
     k = 0
     write (*,'("The MAP estimation of presicion matrix from the data set A")')
     do i=1,ndim
        do j=i,ndim
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
        do j=i,ndim
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
        do j=i,ndim                                                             
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
        do j=i,ndim                                                                
           k = k + 1                                                               
           write (*,'(F8.3,1X)',advance='no'), iAb(i,j)                            
           if (mod(k,5) .eq. 0) then                                               
              write(*,'(/)',advance='no')                                          
           end if                                                                  
        end do                                                                     
     end do                                                                        
     write(*,'(/)',advance='no')                                                   
  end if
    
  write (*,'("The Sparse structure of data set A")')
  do i=1,ndim,1
     do j=i+1,ndim,1
        if ( Aa(i,j) /= 0.0d0 ) then
           write (*,'(I3,1X,I3)'), i,j
        end if
     end do
  end do

  write (*,'("The Sparse structure of data set B")')
  do i=1,ndim,1
     do j=i+1,ndim,1
        if ( Ab(i,j) /= 0.0d0 ) then
           write (*,'(I3,1X,I3)'), i,j
        end if
     end do
  end do

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
    character(20), intent(in) :: programname
    print *,"Usage: ",programname," filename1(data A) filename2(data B)"
    print *,"(option)"
    print *," -mode     :  if 0 normal standarization will be done, if 1 not"
    print *," -ninia    :  initial step of data A"
    print *," -ninib    :  initial step of data B"
    print *," -nfina    :  final   step of data A"
    print *," -nfinb    :  final   step of data B"
    print *," -Sparsity :  sparsity"
    print *," -ndim     :  dimensionality of data"
    
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
       read (unitnum, '(F8.31X)', advance="NO",iostat=err, end=999), v
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
999 N = i !N = i - 1
!  do i = 1,Nndx,1
!     write(*,'(I4,1X)', advance="NO"), atom_list(i)
!     if (mod(i,10) .eq. 0) then
!        write(*,"()")
!     end if
!  end do
  end subroutine read_dat

end program a_ssl
