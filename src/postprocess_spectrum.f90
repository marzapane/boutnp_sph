! =============================================================================
! Computes the instantaneous azimuthal (m) energy spectrum from simulation output.
! A uniform spherical grid is assumed:
!   - Azimuthal  (i=1..n1): n1m=n1-1 unique points, i=n1 is identical to i=1.
!   - Colatitude (k=1..n3): uniform from 0 (N pole) to pi (S pole).
!   - Radial     (j=1..n2): uniform in [r_in, r_in+1]  (gap width = 1).
!
! Usage:
!   ./compute_spectrum <timestep_id> -r=<r_in> [-c|-p]
!
!   -r=<r_in>     : inner radius of the spherical shell (required).
!                   r_ext is deduced as r_in + 1  (gap width = 1).
!   -c (default)  : read Cartesian velocity field from dati/boum*.q
!   -p            : read polar velocity field from dati/spbou*.q and convert
!
! Example:
!   ./compute_spectrum 001000 -r=7.0 -c
!   ./compute_spectrum 001000 -r=7.0 -p
! =============================================================================
! To obtain spectra in physical units multiply by (r_in_real / r_in_code)^3.
! =============================================================================
program compute_spectrum
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  implicit none

  integer :: n1, n2, n3, n1m, n1mh, n3m
  integer :: dummy_int, io_status, arg_count, i, j, k
  real(kind=real64) :: epsil, re, sim_time, pi

  ! Generic velocity arrays (v1,v2,v3) for reading; Cartesian (vx,vy,vz) for computation
  real(kind=real64), allocatable, dimension(:,:,:) :: E_field, v1, v2, v3, pr
  real(kind=real64), allocatable, dimension(:,:,:) :: vx, vy, vz
  real(kind=real64), allocatable, dimension(:) :: r_1d, th_1d
  real(kind=real64), allocatable, dimension(:,:) :: dV, spec, mode
  real(kind=real64) :: dph, dr, dth
  real(kind=real64) :: rint, rext           ! inner and outer radii
  ! Per-point rotation scalars (for polar conversion)
  real(kind=real64) :: colat, lon, cp, sp, ct, st

  character(len=256) :: data_in_file, time_step_str, out_file, coord_type, r_arg
  integer, parameter :: in_unit = 10, out_unit = 11

  ! FFT arrays (fft99 support)
  real(kind=real64), allocatable, dimension(:,:) :: xxr, wwork
  real(kind=real64), allocatable, dimension(:) :: trigxv
  integer, allocatable, dimension(:) :: ifxv
  integer :: id, ip

  ! --------------------------------------------------------------------------
  ! 1. Parse Command Line
  ! --------------------------------------------------------------------------
  arg_count = command_argument_count()
  if (arg_count < 2 .or. arg_count > 3) then
    write(error_unit, '("Usage: ./compute_spectrum <timestep_id> -r=<r_in> [-c|-p]")')
    write(error_unit, '(" -r=<r_in>    : inner radius (required); r_ext = r_in + 1")')
    write(error_unit, '(" -c (default) : Cartesian velocities  -> reads dati/boum*.q")')
    write(error_unit, '(" -p           : Polar velocities       -> reads dati/spbou*.q + converts")')
    write(error_unit, '("Example: ./compute_spectrum 001000 -r=7.0 -p")')
    stop 1
  end if

  call get_command_argument(1, time_step_str)

  ! Parse -r=<value> (second argument, required)
  call get_command_argument(2, r_arg)
  if (r_arg(1:3) /= '-r=') then
    write(error_unit, '("Error: second argument must be -r=<r_in>.")')
    stop 1
  end if
  read(r_arg(4:), *, iostat=io_status) rint
  if (io_status /= 0) then
    write(error_unit, '("Error: could not parse r_in from: ", A)') trim(r_arg)
    stop 1
  end if
  rext = rint + 1.0_real64
  write(*, '(A, F10.4, A, F10.4)') "   Inner radius: ", rint, "   Outer radius: ", rext

  ! Parse optional coordinate flag (third argument)
  if (arg_count == 3) then
    call get_command_argument(3, coord_type)
  else
    coord_type = '-c'
  end if

  ! Choose input file based on coordinate type
  if (trim(coord_type) == '-c') then
    data_in_file = 'dati/boum' // trim(time_step_str) // '.q'
    write(*,*) "-> Coordinate mode: Cartesian"
  else if (trim(coord_type) == '-p') then
    data_in_file = 'dati/spbou' // trim(time_step_str) // '.q'
    write(*,*) "-> Coordinate mode: Polar (will convert to Cartesian)"
  else
    write(error_unit, '("Error: Unknown flag. Must be ''-c'' or ''-p''.")')
    stop 1
  end if


  ! --------------------------------------------------------------------------
  ! 2. Read Data File
  ! --------------------------------------------------------------------------
  write(*,*) "-> Reading data: ", trim(data_in_file)
  open(unit=in_unit, file=trim(data_in_file), form='unformatted', status='old', action='read', iostat=io_status)
  if (io_status /= 0) then
    write(error_unit, '("Error: cannot open data file: ", A)') trim(data_in_file)
    stop 1
  end if

  read(in_unit) dummy_int
  read(in_unit) n1, n2, n3
  write(*, '(a, 3(i0, a))') "   Grid dimensions: (", n1, ", ", n2, ", ", n3, ")"
  read(in_unit) epsil, epsil, re, sim_time

  allocate(E_field(n1, n2, n3), v1(n1, n2, n3), v2(n1, n2, n3), v3(n1, n2, n3), pr(n1, n2, n3))
  ! Read all 5 fields from the single record written by outpf/bounn.f
  read(in_unit) E_field, v1, v2, v3, pr
  close(in_unit)

  n1m = n1 - 1   ! unique azimuthal points (i=n1 is periodic copy of i=1)
  n3m = n3 - 1   ! unique colatitude interior points (poles at k=1 and k=n3)
  n1mh = n1m / 2 + 1
  pi = 2.0_real64 * asin(1.0_real64)

  ! --------------------------------------------------------------------------
  ! 4. Convert to Cartesian if necessary
  ! --------------------------------------------------------------------------
  allocate(vx(n1, n2, n3), vy(n1, n2, n3), vz(n1, n2, n3))

  if (trim(coord_type) == '-c') then
    write(*,*) "-> Using Cartesian velocity components (vx, vy, vz) directly."
    vx = v1; vy = v2; vz = v3

  else ! -p: convert polar (vt, vr, vp) -> Cartesian (vx, vy, vz)
    ! Uniform grid:
    !   lon(i)   = (i-1) * 2*pi / n1m,  i = 1..n1  (i=n1 repeats i=1)
    !   colat(k) = (k-1) * pi  / n3m,   k = 1..n3  (poles at k=1 and k=n3)
    write(*,*) "-> Converting polar velocity components to Cartesian..."
    do k = 1, n3
      colat = real(k-1, real64) * pi / real(n3m, real64)
      cp = cos(colat); sp = sin(colat)
      do j = 1, n2
        do i = 1, n1
          lon = real(i-1, real64) * 2.0_real64 * pi / real(n1m, real64)
          ct = cos(lon); st = sin(lon)
          ! v1=vt (azimuthal), v2=vr (radial), v3=vp (colatitudinal)
          vx(i,j,k) = -v1(i,j,k)*st + v2(i,j,k)*ct*sp + v3(i,j,k)*ct*cp
          vy(i,j,k) =  v1(i,j,k)*ct + v2(i,j,k)*st*sp + v3(i,j,k)*st*cp
          vz(i,j,k) =                  v2(i,j,k)*cp    - v3(i,j,k)*sp
        end do
      end do
    end do

  end if

  ! Free original velocity arrays (no longer needed)
  deallocate(v1, v2, v3, E_field, pr)

  ! --------------------------------------------------------------------------
  ! 5. Build uniform 1D coordinate arrays
  ! --------------------------------------------------------------------------
  ! Colatitude: k=1 -> theta=0 (N pole), k=n3 -> theta=pi (S pole)
  ! Radial: uniform in [r_in, r_ext] with gap width r_ext - r_in = 1
  allocate(r_1d(n2), th_1d(n3))

  do k = 1, n3
    th_1d(k) = real(k-1, real64) * pi / real(n3m, real64)
  end do
  do j = 1, n2
    r_1d(j) = rint + real(j-1, real64) / real(n2-1, real64)   ! uniform in [r_in, r_in+1]
  end do

  ! Uniform grid spacings
  dr  = 1.0_real64 / real(n2-1, real64)   ! dr = (r_ext - r_in) / (n2-1) = 1/(n2-1)
  dth = pi / real(n3m, real64)             ! d(theta) = pi/(n3-1)
  dph = 2.0_real64 * pi                   ! full azimuthal integral

  ! Volume element matrix dV(j, k) = r^2 * sin(theta) * dr * d(theta) * d(phi)
  ! Note: dV = 0 at poles (sin(0) = sin(pi) = 0), so pole points contribute nothing.
  allocate(dV(n2, n3))
  do k = 1, n3
    do j = 1, n2
      dV(j, k) = r_1d(j)**2 * sin(th_1d(k)) * dr * dth * dph
    end do
  end do

  ! --------------------------------------------------------------------------
  ! 6. Compute Azimuthal (m) Spectrum using fft99
  ! --------------------------------------------------------------------------
  write(*,*) "-> Computing Energy Spectrum using fft99..."

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

  ! Loop over radial shells; accumulate energy from all three Cartesian components
  do j = 1, n2

    ! --- vx ---
    xxr = 0.0_real64
    do k = 1, n3
      xxr(1, k) = vx(n1m, j, k)
      do i = 1, n1m
        xxr(i+1, k) = vx(i, j, k)
      end do
      xxr(n1m+2, k) = vx(1, j, k)
    end do
    call fft99(xxr, wwork, trigxv, ifxv, 1, n1m+2, n1m, n3, -1)
    do i = 1, n1mh
      ip = 2*i
      id = 2*i - 1
      do k = 1, n3
        spec(k, i) = spec(k, i) + 0.5_real64 * (xxr(id,k)**2 + xxr(ip,k)**2) * dV(j, k)
      end do
    end do

    ! --- vy ---
    xxr = 0.0_real64
    do k = 1, n3
      xxr(1, k) = vy(n1m, j, k)
      do i = 1, n1m
        xxr(i+1, k) = vy(i, j, k)
      end do
      xxr(n1m+2, k) = vy(1, j, k)
    end do
    call fft99(xxr, wwork, trigxv, ifxv, 1, n1m+2, n1m, n3, -1)
    do i = 1, n1mh
      ip = 2*i
      id = 2*i - 1
      do k = 1, n3
        spec(k, i) = spec(k, i) + 0.5_real64 * (xxr(id,k)**2 + xxr(ip,k)**2) * dV(j, k)
      end do
    end do

    ! --- vz ---
    xxr = 0.0_real64
    do k = 1, n3
      xxr(1, k) = vz(n1m, j, k)
      do i = 1, n1m
        xxr(i+1, k) = vz(i, j, k)
      end do
      xxr(n1m+2, k) = vz(1, j, k)
    end do
    call fft99(xxr, wwork, trigxv, ifxv, 1, n1m+2, n1m, n3, -1)
    do i = 1, n1mh
      ip = 2*i
      id = 2*i - 1
      do k = 1, n3
        spec(k, i) = spec(k, i) + 0.5_real64 * (xxr(id,k)**2 + xxr(ip,k)**2) * dV(j, k)
      end do
    end do

  end do

  ! Apply factor of 2 for negative frequencies (conjugate folding), except m=0 and Nyquist
  do i = 2, n1mh
    if (i == n1mh .and. mod(n1m, 2) == 0) cycle  ! Nyquist: no factor 2
    spec(:, i) = spec(:, i) * 2.0_real64
  end do

  ! --------------------------------------------------------------------------
  ! 7. Output the spectrum
  ! --------------------------------------------------------------------------
  out_file = 'spctr' // trim(time_step_str) // '.bin'
  write(*,*) "-> Writing output to ", trim(out_file)
  open(unit=out_unit, file=trim(out_file), form='unformatted', status='replace')
  write(out_unit) n3m, n1mh
  write(out_unit) ((mode(k, i), k=1, n3m), i=1, n1mh)
  write(out_unit) ((spec(k, i), k=1, n3m), i=1, n1mh)
  write(out_unit) sim_time
  close(out_unit)

  write(*,*) "-> Done."

end program compute_spectrum
