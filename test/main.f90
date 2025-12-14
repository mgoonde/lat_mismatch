program main
  use routines
  implicit none

  integer :: ios
  character(len=128) :: msg
  real(rp) :: a0_a, a0_b
  real(rp) :: a1_a(3), a2_a(3), a3_a(3)
  real(rp) :: a1_b(3), a2_b(3), a3_b(3)
  integer :: nat_a, nat_b
  real(rp), allocatable :: bas_a(:,:), bas_b(:,:)
  real(rp) :: lat_a(3,3), lat_b(3,3), invlat_a(3,3), invlat_b(3,3)
  integer :: h, k, ierr

  namelist/unitcell/a0_a, a0_b, a1_a, a1_b, a2_a, a2_b, a3_a, a3_b, nat_a, nat_b
  namelist/basis/bas_a, bas_b

  ! read unitcell
  read(*, nml=unitcell, iostat=ios, iomsg=msg )
  if( ios /= 0 ) then
     write(*,*) "err reading nml=unitcell, ios=",ios
     write(*,*) "msg:", trim(msg)
     error stop 1
  end if
  ! alloc
  allocate( bas_a(1:3,1:nat_a) )
  allocate( bas_b(1:3,1:nat_b) )
  ! read basis vecs
  read(*, nml=basis, iostat=ios, iomsg=msg )
  if( ios /= 0 ) then
     write(*,*) "err reading nml=basis, ios=",ios
     write(*,*) "msg:", trim(msg)
     error stop 1
  end if


  ! lat in columns
  lat_a(:,1) = a1_a
  lat_a(:,2) = a2_a
  lat_a(:,3) = a3_a
  lat_b(:,1) = a1_b
  lat_b(:,2) = a2_b
  lat_b(:,3) = a3_b

  call invmat3x3( lat_a, invlat_a )
  call invmat3x3( lat_b, invlat_b )

  write(*,*) "sys A:"
  write(*,*) "a0:",a0_a
  write(*,*) "lat vecs:"
  write(*,*) lat_a(:,1)
  write(*,*) lat_a(:,2)
  write(*,*) lat_a(:,3)
  write(*,*) "bas vecs:"
  write(*,"(3f9.5)") bas_a

  write(*,*)
  write(*,*) "sys B:"
  write(*,*) "a0:",a0_b
  write(*,*) "lat vecs:"
  write(*,*) lat_b(:,1)
  write(*,*) lat_b(:,2)
  write(*,*) lat_b(:,3)
  write(*,*) "bas vecs:"
  write(*,"(3f9.5)") bas_b

  ierr = frac_approx( a0_a/a0_b, h, k, kmax=50, nmin=3 )
  if( ierr /= 0 ) error stop 1

  write(*,"(a,1x,f9.6,1x,a,1x,i0,'/',i0)") "frac:",a0_a/a0_b,"has approx integer fractional h/k:",h,k

  write(*,*) " repetitions:"
  write(*,"(1x,'k*a0_b = ',i0,'*',g0.6,' = ',f12.6)") k, a0_a, k*a0_a
  write(*,"(1x,'k*a0_b = ',i0,'*',g0.6,' = ',f12.6)") h, a0_b, h*a0_b
  write(*,*) "error = ", (k*a0_a - h*a0_b)/(k*a0_a)

end program main
