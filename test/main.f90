program main
  use routines
  use, intrinsic :: iso_fortran_env, only: io_end => iostat_end
  implicit none

  integer :: ios, u0
  character(len=128) :: msg
  real(rp) :: a0_a, a0_b
  real(rp) :: a1_a(3), a2_a(3), a3_a(3)
  real(rp) :: a1_b(3), a2_b(3), a3_b(3)
  integer :: nat_a, nat_b
  real(rp), allocatable :: bas_a(:,:), bas_b(:,:)
  real(rp) :: lat_a(3,3), lat_b(3,3), invlat_a(3,3), invlat_b(3,3)
  integer, allocatable :: typ_a(:), typ_b(:)
  integer :: h, k, ierr
  integer :: kmax, nmin
  integer :: xx, yy, zz, i
  integer :: n_layers_a, n_layers_b
  real(rp) :: interlattice_dist
  real(rp) :: shift(3), dz(3), vacuum
  character(*), parameter :: fout="tt.xyz"

  namelist/unitcell/a0_a, a0_b, a1_a, a1_b, a2_a, a2_b, a3_a, a3_b, nat_a, nat_b
  namelist/basis/bas_a, bas_b, typ_a, typ_b
  namelist/frac/kmax, nmin
  namelist/z_spacing/n_layers_a, n_layers_b, interlattice_dist, vacuum

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
  allocate( typ_a(1:nat_a), source=1)
  allocate( typ_b(1:nat_b), source=1)
  ! read basis vecs
  read(*, nml=basis, iostat=ios, iomsg=msg )
  if( ios /= 0 ) then
     write(*,*) "err reading nml=basis, ios=",ios
     write(*,*) "msg:", trim(msg)
     error stop 2
  end if

  ! read the frac nml
  kmax = 50
  nmin = 3
  read(*, nml=frac, iostat=ios, iomsg=msg )
  if( ios /= 0 .and. ios /= io_end ) then
     write(*,*) "err reading nml=frac, ios=",ios
     write(*,*) "msg:",trim(msg)
     error stop 3
  end if

  ! read the z_spacing nml
  interlattice_dist = 0.0_rp
  vacuum = 10.0_rp
  n_layers_a = 1; n_layers_b = 1
  read(*, nml=z_spacing, iostat=ios, iomsg=msg )
  if( ios /= 0 ) then
     write(*,*) "err reading nml=z_spacing, ios=",ios
     write(*,*) "msg:", trim(msg)
     error stop 4
  end if

  ! unitcell lat in columns
  lat_a(:,1) = a1_a
  lat_a(:,2) = a2_a
  lat_a(:,3) = a3_a
  lat_a = lat_a*a0_a
  lat_b(:,1) = a1_b
  lat_b(:,2) = a2_b
  lat_b(:,3) = a3_b
  lat_b = lat_b*a0_b

  ! call invmat3x3( lat_a, invlat_a )
  ! call invmat3x3( lat_b, invlat_b )

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

  ierr = frac_approx( a0_a/a0_b, h, k, kmax=kmax, nmin=nmin )
  if( ierr /= 0 ) error stop 1

  write(*,*)
  write(*,"(a,1x,f9.6,1x,a,1x,i0,'/',i0)") "frac:",a0_a/a0_b,"has approx integer fractional h/k:",h,k

  write(*,*) " repetitions:"
  write(*,"(1x,'k*a0_b = ',i0,'*',g0.6,' = ',f12.6)") k, a0_a, k*a0_a
  write(*,"(1x,'k*a0_b = ',i0,'*',g0.6,' = ',f12.6)") h, a0_b, h*a0_b
  write(*,"(1x,a,1x,e0.6)") "absolute error =", abs(k*a0_a - h*a0_b)
  write(*,"(1x,a,1x,e0.6)") "relative error =", abs( (k*a0_a - h*a0_b)/(k*a0_a) )


  ! shift for system b, such that it starts at same value of z as last layer of a
  dz = [0.0_rp, 0.0_rp, real(n_layers_a)-0.5_rp]
  dz = matmul(lat_a, dz)

  open(newunit=u0, file=fout, status="replace")
  write(u0, *) nat_a*k*k*n_layers_a + nat_b*h*h*n_layers_b
  write(u0,*) 'Lattice="', lat_a(:,1)*k, lat_a(:,2)*k, &
       lat_a(:,3)*(n_layers_a) &
       + lat_b(:,3)*n_layers_b &
       - matmul(lat_a,[0.0_rp, 0.0_rp, 0.5_rp]) &
       + [0.0_rp, 0.0_rp, interlattice_dist + vacuum] &
       , '"'
  ! write sys a, which is repeated k times in xy, and n_layers_a in z
  do xx = 0, k-1
     do yy = 0, k-1
        do zz = 0, n_layers_a-1
           shift = real( [xx, yy, zz])
           do i = 1, nat_a
              write(u0,*) typ_a(i), matmul(lat_a, shift+bas_a(:,i))
           end do
        end do
     end do
  end do
  ! write sys b, which is repeated h times in xy, and n_layers_a:n_layers_b in z
  do xx = 0, h-1
     do yy = 0, h-1
        do zz = 0, n_layers_b-1
           shift = real([xx,yy,zz])
           do i = 1, nat_b
                   ! + matmul(lat_b, [0.0_rp, 0.0_rp, interlayer_dist])
              write(u0,*) typ_b(i), matmul(lat_b, shift+bas_b(:,i)) + dz + [0.0_rp, 0.0_rp, interlattice_dist]
           end do
        end do
     end do
  end do
  close(u0, status="keep")

  write(*,*) " >> Output written to: "//trim(fout)

end program main
