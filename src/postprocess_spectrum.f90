! =============================================================================
! Computes the instantaneous surface spectrum (radially averaged) from output radial field in spbou*.q
! =============================================================================
program compute_spectrum
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  implicit none

  integer :: n1, n2, n3, n1m, n1mh, n1g, n2g, n3g
  integer :: dummy_int, io_status, arg_count, i, j, k
  real(kind=real64) :: epsil, re, sim_time, pi

  real(kind=real64), allocatable, dimension(:,:,:) :: E_field, vt, vr, vp, pr
  real(kind=real64), allocatable, dimension(:,:,:) :: xxx, yyy, zzz
  real(kind=real64), allocatable, dimension(:) :: r_1d, th_1d, dr_1d, dth_1d
  real(kind=real64), allocatable, dimension(:,:) :: dV, spec, mode
  real(kind=real64) :: dph

  character(len=256) :: data_in_file, grid_in_file, time_step_str, out_file
  integer, parameter :: in_unit = 10, out_unit = 11

  ! FFT arrays (fft99 support)
  real(kind=real64), allocatable, dimension(:,:) :: xxr, wwork
  real(kind=real64), allocatable, dimension(:) :: trigxv
  integer, allocatable, dimension(:) :: ifxv
  integer :: id, ip

  ! Check arguments
  arg_count = command_argument_count()
  if (arg_count /= 1) then
    write(error_unit, '("Usage: ./compute_spectrum <timestep_id>")')
    write(error_unit, '("Example: ./compute_spectrum 000100")')
    stop 1
  end if

  call get_command_argument(1, time_step_str)
  data_in_file = 'dati/spbou' // trim(time_step_str) // '.q'
  grid_in_file = 'grid.xyz'

  ! Read data file
  write(*,*) "-> Reading ", trim(data_in_file)
  open(unit=in_unit, file=trim(data_in_file), form='unformatted', status='old', action='read', iostat=io_status)
  if (io_status /= 0) stop "Error: cannot open data input file."
  
  read(in_unit) dummy_int
  read(in_unit) n1, n2, n3
  write(*, '(a, 3(i0, a))') "   Grid dimensions: (", n1, ", ", n2, ", ", n3, ")"
  read(in_unit) epsil, epsil, re, sim_time
  
  allocate(E_field(n1, n2, n3), vt(n1, n2, n3), vr(n1, n2, n3), vp(n1, n2, n3), pr(n1, n2, n3))
  read(in_unit) E_field, vt, vr, vp, pr
  close(in_unit)

  ! Read grid
  write(*,*) "-> Reading ", trim(grid_in_file)
  open(unit=in_unit, file=trim(grid_in_file), form='unformatted', status='old', action='read', iostat=io_status)
  if (io_status /= 0) stop "Error: cannot open grid file."

  read(in_unit) dummy_int
  read(in_unit) n1g, n2g, n3g
  if (n1g /= n1 .or. n2g /= n2 .or. n3g /= n3) then
    stop "Error: grid and data file dimensions do not correspond."
  end if
  allocate(xxx(n1, n2, n3), yyy(n1, n2, n3), zzz(n1, n2, n3))
  read(in_unit) xxx, yyy, zzz
  close(in_unit)

  n1m = n1 - 1
  n1mh = n1m / 2 + 1
  pi = 2.0_real64 * asin(1.0_real64)

  ! -------------------------------------------------------------
  ! Compute 1D profiles for coordinates r and phi by averaging
  ! -------------------------------------------------------------
  allocate(r_1d(n2), th_1d(n3))
  
  r_1d = 0.0_real64
  do k = 1, n3
    do j = 1, n2
      do i = 1, n1
        r_1d(j) = r_1d(j) + sqrt(xxx(i,j,k)**2 + yyy(i,j,k)**2 + zzz(i,j,k)**2)
      end do
    end do
  end do
  r_1d = r_1d / real(n1 * n3, real64)

  th_1d = 0.0_real64
  do k = 1, n3
    do j = 1, n2
      do i = 1, n1
        th_1d(k) = th_1d(k) + acos(zzz(i,j,k) / sqrt(xxx(i,j,k)**2 + yyy(i,j,k)**2 + zzz(i,j,k)**2))
      end do
    end do
    th_1d(k) = th_1d(k) / real(n1 * n2, real64)
  end do

  ! -------------------------------------------------------------
  ! Compute differentials using central differences
  ! -------------------------------------------------------------
  allocate(dr_1d(n2), dth_1d(n3))
  dr_1d = 0.0_real64
  dth_1d = 0.0_real64
  
  if (n2 > 1) then
    dr_1d(1) = r_1d(2) - r_1d(1)
    do j = 2, n2-1
      dr_1d(j) = (r_1d(j+1) - r_1d(j-1)) / 2.0_real64
    end do
    dr_1d(n2) = r_1d(n2) - r_1d(n2-1)
  end if

  if (n3 > 1) then
    dth_1d(1) = th_1d(2) - th_1d(1)
    do k = 2, n3-1
      dth_1d(k) = (th_1d(k+1) - th_1d(k-1)) / 2.0_real64
    end do
    dth_1d(n3) = th_1d(n3) - th_1d(n3-1)
  end if

  dph = 2.0_real64 * pi / real(n1m, real64)

  ! Volume element matrix dV(j, k)
  allocate(dV(n2, n3))
  do k = 1, n3
    do j = 1, n2
      dV(j, k) = (r_1d(j)**2) * sin(th_1d(k)) * dr_1d(j) * dth_1d(k) * dph
    end do
  end do

  ! -------------------------------------------------------------
  ! Compute Spectrum using fft99
  ! -------------------------------------------------------------
  write(*,*) "-> Computing Energy Spectrum using fft99..."
  
  ! Setup fftfax
  allocate(trigxv(3*n1m/2 + 1), ifxv(13))
  allocate(xxr(n1m+2, n3), wwork(n1m+1, n3))
  call fftfax(n1m, ifxv, trigxv)

  allocate(spec(n3, n1mh), mode(n3, n1mh))
  spec = 0.0_real64
  mode = 0.0_real64

  do i = 1, n1mh
    do k = 1, n3
      mode(k, i) = real(i-1, real64) / max(sin(th_1d(k)), 1e-12_real64)
    end do
  end do

  ! Loop over radial points to compute volume integral of energy
  do j = 1, n2
    
    ! --- Azimuthal velocity (vt) ---
    xxr = 0.0_real64
    do k = 1, n3
      xxr(1, k) = vt(n1m, j, k)
      do i = 1, n1m
        xxr(i+1, k) = vt(i, j, k)
      end do
      xxr(n1m+2, k) = vt(1, j, k)
    end do
    ! Perform parallel FFT along n3 directions
    call fft99(xxr, wwork, trigxv, ifxv, 1, n1m+2, n1m, n3, -1)
    
    do i = 1, n1mh
      ip = 2*i
      id = 2*i - 1
      do k = 1, n3
        spec(k, i) = spec(k, i) + 0.5_real64 * (xxr(id, k)**2 + xxr(ip, k)**2) * dV(j, k)
      end do
    end do

    ! --- Radial velocity (vr) ---
    xxr = 0.0_real64
    do k = 1, n3
      xxr(1, k) = vr(n1m, j, k)
      do i = 1, n1m
        xxr(i+1, k) = vr(i, j, k)
      end do
      xxr(n1m+2, k) = vr(1, j, k)
    end do
    call fft99(xxr, wwork, trigxv, ifxv, 1, n1m+2, n1m, n3, -1)
    
    do i = 1, n1mh
      ip = 2*i
      id = 2*i - 1
      do k = 1, n3
        spec(k, i) = spec(k, i) + 0.5_real64 * (xxr(id, k)**2 + xxr(ip, k)**2) * dV(j, k)
      end do
    end do

    ! --- Polar velocity (vp) ---
    xxr = 0.0_real64
    do k = 1, n3
      xxr(1, k) = vp(n1m, j, k)
      do i = 1, n1m
        xxr(i+1, k) = vp(i, j, k)
      end do
      xxr(n1m+2, k) = vp(1, j, k)
    end do
    call fft99(xxr, wwork, trigxv, ifxv, 1, n1m+2, n1m, n3, -1)
    
    do i = 1, n1mh
      ip = 2*i
      id = 2*i - 1
      do k = 1, n3
        spec(k, i) = spec(k, i) + 0.5_real64 * (xxr(id, k)**2 + xxr(ip, k)**2) * dV(j, k)
      end do
    end do
  end do

  ! Apply negative frequency factor (factor of 2 for m > 0, except Nyquist if n1m is even)
  ! We apply this factor since fft99 yields C(k) coefficients that sum over full range [-N/2, N/2]
  ! and we are folding negative frequencies into positive ones
  do i = 2, n1mh
    if (i == n1mh .and. mod(n1m, 2) == 0) then
      ! Nyquist mode, no factor 2
      continue
    else
      spec(:, i) = spec(:, i) * 2.0_real64
    end if
  end do

  ! -------------------------------------------------------------
  ! Output the spectrum
  ! -------------------------------------------------------------
  out_file = 'spctr' // trim(time_step_str) // '.bin'
  write(*,*) "-> Writing output to ", trim(out_file)
  open(unit=out_unit, file=trim(out_file), form='unformatted', status='replace')
  write(out_unit) n3, n1mh
  write(out_unit) ((mode(k, i), k=1, n3), i=1, n1mh)
  write(out_unit) ((spec(k, i), k=1, n3), i=1, n1mh)
  write(out_unit) sim_time
  close(out_unit)

  write(*,*) "-> Done."

end program compute_spectrum
