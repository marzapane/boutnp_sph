! =============================================================================
! SHT_3D_POSTPROCESS
! 
! This utility performs a Spherical Harmonic Transform (SHT) on a 3D vector 
! velocity field using the SHTOOLS library.
!
! Features:
! - Modular design separating I/O and SHT processing.
! - Computes the total kinetic energy spectrum, and separate m=0 / m!=0 spectra.
! - Properly handles both polar and Cartesian input velocity fields.
!
! Compilation (example):
! gfortran -O3 sht_3d_postprocess.f90 -o sht_3d_utility -lshtools -lfftw3 -lm
! gfortran -O3 src/sht_3d_postprocess.f90 -o sht_3d_utility -I/usr/local/include -m64 -fPIC -O3 -std=gnu -ffast-math -L/usr/local/lib -lSHTOOLS -L/usr/local -lfftw3 -lm -llapack -lblas
! =============================================================================

module data_io_mod
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

contains

    subroutine read_simulation_data(file_path, n1, n2, n3, sim_time, v1, v2, v3)
        character(len=*), intent(in) :: file_path
        integer, intent(out) :: n1, n2, n3
        real(kind=real64), intent(out) :: sim_time
        real(kind=real64), allocatable, intent(out) :: v1(:,:,:), v2(:,:,:), v3(:,:,:)
        
        integer :: in_unit = 10, io_status, dummy_int
        real(kind=real64) :: epsil, re
        real(kind=real64), allocatable :: field_temp(:,:,:), pr(:,:,:)
        
        write(*,*) "-> Reading data from: ", trim(file_path)
        open(unit=in_unit, file=trim(file_path), form='unformatted', status='old', action='read', iostat=io_status)
        if (io_status /= 0) then
            write(error_unit, '("Error: Could not open ", A)') trim(file_path)
            stop 1
        end if

        read(in_unit) dummy_int
        read(in_unit) n1, n2, n3
        read(in_unit) epsil, epsil, re, sim_time
        
        allocate(field_temp(n1, n2, n3), pr(n1, n2, n3))
        allocate(v1(n1, n2, n3), v2(n1, n2, n3), v3(n1, n2, n3))
        
        ! Read all fields (standard output format from outpf)
        read(in_unit) field_temp ! Velocity Magnitude
        read(in_unit) v1         ! Azimuthal (vt) or X (vx)
        read(in_unit) v2         ! Radial (vr) or Y (vy)
        read(in_unit) v3         ! Polar (vp) or Z (vz)
        read(in_unit) pr         ! Pressure
        close(in_unit)
        
    end subroutine read_simulation_data

end module data_io_mod


module sht_processing_mod
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    use SHTOOLS, only: SHExpandDH, SHPowerSpectrum
    implicit none

contains

    ! 1. Computes the transform and returns the spherical harmonics coefficients
    subroutine get_sht_coeffs(grid, n_lat, lmax_calc, coeffs)
        integer, intent(in) :: n_lat, lmax_calc
        integer             :: lmax
        real(kind=real64), intent(in) :: grid(n_lat, 2*n_lat)
        real(kind=real64), intent(out) :: coeffs(2, lmax_calc+1, lmax_calc+1)
        
        ! Perform SHT (Grid -> Coefficients)
        ! sampling=2: grid is equally spaced (n x 2n)
        call SHExpandDH(grid, n_lat, coeffs, lmax, lmax_calc=lmax_calc, sampling=2)
    end subroutine get_sht_coeffs

    ! 2. Computes and returns the energy spectrum
    subroutine get_energy_spectrum(coeffs, lmax, spec_tot)
        integer, intent(in) :: lmax
        real(kind=real64), intent(in) :: coeffs(2, lmax+1, lmax+1)
        real(kind=real64), intent(out) :: spec_tot(lmax+1)
        
        ! Compute Power Spectrum
        call SHPowerSpectrum(coeffs, lmax, spec_tot)
    end subroutine get_energy_spectrum

    ! 3. Computes two separate spectra, one for the m=0 mode and one for all the rest
    subroutine get_m_spectra(coeffs, lmax, spec_m0, spec_m_rest)
        integer, intent(in) :: lmax
        real(kind=real64), intent(in) :: coeffs(2, lmax+1, lmax+1)
        real(kind=real64), intent(out) :: spec_m0(lmax+1), spec_m_rest(lmax+1)
        
        real(kind=real64) :: spec_tot(lmax+1)
        integer :: l
        
        ! Get the total power spectrum first
        call SHPowerSpectrum(coeffs, lmax, spec_tot)
        
        do l = 0, lmax
            ! The power at m=0 is just the square of the (1, l+1, 1) coefficient
            ! (assuming default 4-pi normalized harmonics)
            spec_m0(l+1) = coeffs(1, l+1, 1)**2
            
            ! The rest is the total minus the m=0 component
            spec_m_rest(l+1) = spec_tot(l+1) - spec_m0(l+1)
        end do
    end subroutine get_m_spectra

end module sht_processing_mod


program sht_3d_utility
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    use data_io_mod
    use sht_processing_mod
    implicit none

    integer :: n1, n2, n3, lmax, n_lat, n_lon
    integer :: arg_count, i, j, k
    real(kind=real64) :: sim_time, pi, lon, colat, ct, st, cp, sp
    
    character(len=256) :: input_path, coord_type, out_file
    integer, parameter :: out_unit = 11
    
    real(kind=real64), allocatable :: v1(:,:,:), v2(:,:,:), v3(:,:,:)
    real(kind=real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
    
    real(kind=real64), allocatable :: grid_x(:,:), grid_y(:,:), grid_z(:,:)
    real(kind=real64), allocatable :: coeffs_x(:,:,:), coeffs_y(:,:,:), coeffs_z(:,:,:)
    real(kind=real64), allocatable :: spec_tot(:), spec_m0(:), spec_m_rest(:)
    real(kind=real64), allocatable :: spec_tot_3d(:,:), spec_m0_3d(:,:), spec_m_rest_3d(:,:)
    
    real(kind=real64) :: vt, vr, vp

    ! 1. Parse Command Line
    arg_count = command_argument_count()
    if (arg_count < 1) then
        write(error_unit, '("Usage: ./sht_3d_utility <file_path> [-c|-p]")')
        write(error_unit, '(" -c (default) = velocity field input is in Cartesian coordinates")')
        write(error_unit, '(" -p           = velocity field input is in polar coordinates")')
        stop 1
    end if

    call get_command_argument(1, input_path)
    
    if (arg_count >= 2) then
        call get_command_argument(2, coord_type)
    else
        coord_type = '-c'
    end if

    ! 2. Read Data
    call read_simulation_data(input_path, n1, n2, n3, sim_time, v1, v2, v3)
    write(*, '(A, 3(I0, A), F10.4)') "   Dimensions: ", n1, "x", n2, "x", n3, " at time ", sim_time

    ! 3. Handle Velocity Field Appropriately (Polar vs Cartesian)
    allocate(vx(n1, n2, n3), vy(n1, n2, n3), vz(n1, n2, n3))
    
    if (trim(coord_type) == '-p') then
        write(*,*) "-> Converting polar velocity components to Cartesian..."
        pi = 2.0_real64 * asin(1.0_real64)
        
        do k = 1, n3
            ! Assuming colatitude goes from 0 to pi uniformly
            colat = real(k-1, real64) * pi / real(n3-1, real64)
            cp = cos(colat)
            sp = sin(colat)
            do j = 1, n2
                do i = 1, n1
                    ! Assuming longitude goes from 0 to 2*pi uniformly
                    lon = real(i-1, real64) * 2.0_real64 * pi / real(n1-1, real64)
                    ct = cos(lon)
                    st = sin(lon)
                    
                    vt = v1(i, j, k)
                    vr = v2(i, j, k)
                    vp = v3(i, j, k)
                    
                    ! Standard rotation matrix for spherical to cartesian vectors
                    vx(i, j, k) = -vt*st + vr*ct*sp + vp*ct*cp
                    vy(i, j, k) =  vt*ct + vr*st*sp + vp*st*cp
                    vz(i, j, k) =          vr*cp    - vp*sp
                end do
            end do
        end do
    else if (trim(coord_type) == '-c') then
        write(*,*) "-> Using Cartesian velocity components directly."
        vx = v1; vy = v2; vz = v3
    else
        write(error_unit, '("Error: Unknown coordinate type option. Must be ''-c'' or ''-p''.")')
        stop 1
    end if
    
    ! We don't need the original arrays anymore
    deallocate(v1, v2, v3)

    ! 4. Setup SHT Parameters
    n_lat = n3 - 1
    n_lon = n1 - 1
    lmax = n_lat / 2 - 1

    allocate(grid_x(n_lat, n_lon), grid_y(n_lat, n_lon), grid_z(n_lat, n_lon))
    allocate(coeffs_x(2, lmax+1, lmax+1), coeffs_y(2, lmax+1, lmax+1), coeffs_z(2, lmax+1, lmax+1))
    allocate(spec_tot(lmax+1), spec_m0(lmax+1), spec_m_rest(lmax+1))
    
    allocate(spec_tot_3d(lmax+1, n2), spec_m0_3d(lmax+1, n2), spec_m_rest_3d(lmax+1, n2))
    
    spec_tot_3d = 0.0_real64
    spec_m0_3d = 0.0_real64
    spec_m_rest_3d = 0.0_real64

    write(*,*) "-> Performing SHT on Velocity Field..."
    do j = 1, n2
        ! Extract 2D spherical slices, transpose to (lat, lon) format for SHTOOLS
        grid_x(:, :) = transpose(vx(1:n_lon, j, 1:n_lat))
        grid_y(:, :) = transpose(vy(1:n_lon, j, 1:n_lat))
        grid_z(:, :) = transpose(vz(1:n_lon, j, 1:n_lat))
        
        ! Call Subroutine 1: Get Coefficients
        call get_sht_coeffs(grid_x, n_lat, lmax, coeffs_x)
        call get_sht_coeffs(grid_y, n_lat, lmax, coeffs_y)
        call get_sht_coeffs(grid_z, n_lat, lmax, coeffs_z)

        ! Call Subroutine 2: Get Energy Spectra
        call get_energy_spectrum(coeffs_x, lmax, spec_tot)
        spec_tot_3d(:, j) = spec_tot_3d(:, j) + spec_tot(:)
        
        call get_energy_spectrum(coeffs_y, lmax, spec_tot)
        spec_tot_3d(:, j) = spec_tot_3d(:, j) + spec_tot(:)
        
        call get_energy_spectrum(coeffs_z, lmax, spec_tot)
        spec_tot_3d(:, j) = spec_tot_3d(:, j) + spec_tot(:)

        ! Call Subroutine 3: Get m=0 and m!=0 Spectra
        call get_m_spectra(coeffs_x, lmax, spec_m0, spec_m_rest)
        spec_m0_3d(:, j) = spec_m0_3d(:, j) + spec_m0(:)
        spec_m_rest_3d(:, j) = spec_m_rest_3d(:, j) + spec_m_rest(:)
        
        call get_m_spectra(coeffs_y, lmax, spec_m0, spec_m_rest)
        spec_m0_3d(:, j) = spec_m0_3d(:, j) + spec_m0(:)
        spec_m_rest_3d(:, j) = spec_m_rest_3d(:, j) + spec_m_rest(:)
        
        call get_m_spectra(coeffs_z, lmax, spec_m0, spec_m_rest)
        spec_m0_3d(:, j) = spec_m0_3d(:, j) + spec_m0(:)
        spec_m_rest_3d(:, j) = spec_m_rest_3d(:, j) + spec_m_rest(:)

        if (mod(j, 25) == 0) write(*,*) "      Progress: shell ", j, "/", n2
    end do

    ! 5. Output the Three Spectra
    out_file = 'spectrum.dat'
    
    write(*,*) "-> Saving Spectrum to: ", trim(out_file)
    open(unit=out_unit, file=trim(out_file), status='replace')
    write(out_unit, '(A)') "# L-degree | Radial Index | Total Energy Spectrum | Zonal Energy Spectrum (m=0) | Residual Energy Spectrum (m!=0)"
    do j = 1, n2
        do i = 0, lmax
            write(out_unit, '(I6, I6, E15.7, E15.7, E15.7)') i, j, spec_tot_3d(i+1, j), spec_m0_3d(i+1, j), spec_m_rest_3d(i+1, j)
        end do
        write(out_unit, *)
    end do
    close(out_unit)

end program sht_3d_utility
