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
        read(in_unit) field_temp, v1, v2, v3, pr
        close(in_unit)

    end subroutine read_simulation_data

end module data_io_mod


module sht_processing_mod
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    use SHTOOLS, only: SHExpandDH, SHPowerSpectrum
    use fitpack, only: fitpack_curve, FITPACK_SUCCESS, FITPACK_MESSAGE, FP_FLAG
    implicit none

    real(kind=real64), parameter, private :: pi = 2.0_real64 * asin(1.0_real64)
    real(kind=real64) :: s_factor_mult = 0.005_real64

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

    ! 2. Computes the energy spectrum, the m=0 spectrum and all the rest,
    ! and accumulates the results.
    subroutine upd_energy_spectrum(coeffs, lmax, spec_tot, spec_m0, spec_rest)
        integer, intent(in) :: lmax
        real(kind=real64), intent(in) :: coeffs(2, 0:lmax, 0:lmax)
        real(kind=real64), intent(inout) :: spec_tot(0:lmax)
        real(kind=real64), intent(inout) :: spec_m0(0:lmax), spec_rest(0:lmax)

        real(kind=real64) :: spec(0:lmax)
        integer :: l

        ! Get the total power spectrum first
        call SHPowerSpectrum(coeffs, lmax, spec)

        spec_tot(:) = spec_tot(:) + spec(:)

        do l = 0, lmax
            ! The power at m=0 is just the square of the (1, l, 0) coefficient
            ! (assuming default 4-pi normalized harmonics)
            spec_m0(l) = spec_m0(l) + coeffs(1, l, 0)**2

            ! The rest is the total minus the m=0 component
            spec_rest(l) = spec_rest(l) + (spec(l) - coeffs(1, l, 0)**2)
        end do
    end subroutine upd_energy_spectrum

    ! 3. Bins a 2D horizontal spectrum (spec(0:lmax, n_shells)) logarithmically
    !    into n_bins. Employs CIC interpolation.
    !    Includes physical volume integration weights (0.5 * spec * r_j^2 * dr).
    subroutine bin_horizontal_spectrum(spec, lmax, n_shells, r_min, r_max, n_bins, E_K, K_bins)
        integer,           intent(in)  :: lmax, n_shells, n_bins
        real(kind=real64), intent(in)  :: spec(0:lmax, n_shells)
        real(kind=real64), intent(in)  :: r_min, r_max
        real(kind=real64), intent(out) :: E_K(0:n_bins), K_bins(0:n_bins)

        real(kind=real64) :: dr, log_K_min, log_K_max, inv_delta_log, delta_log
        real(kind=real64) :: r_j, K_tot, E_contrib, real_idx, w_lo, w_hi
        integer :: l, j, idx

        dr = (r_max - r_min) / real(n_shells - 1, real64)

        ! Use identical bin boundaries as the 3D volume spectrum for direct comparability
        log_K_min     = log(min(sqrt(2.0_real64)/r_max, 2*pi/(r_max-r_min)))
        log_K_max     = log(real(lmax*(lmax+1), real64)/r_min**2 &
                  + (2*sqrt(real(n_shells*(n_shells+1), real64))*2*pi/(r_max-r_min))**2) / 2.0_real64
        delta_log     = (log_K_max - log_K_min) / real(n_bins, real64)
        inv_delta_log = real(n_bins, real64) / (log_K_max - log_K_min)

        ! Build logarithmic bin centres
        K_bins(0) = 0.0_real64
        do j = 1, n_bins
            K_bins(j) = exp(log_K_min + (real(j, real64) - 0.5_real64) * delta_log)
        end do

        E_K = 0.0_real64

        ! ---- CIC accumulation ----
        do j = 1, n_shells
            r_j = r_min + real(j - 1, real64) * dr
            do l = 0, lmax
                ! Physical energy density: 1/2 * (4*pi) * spec * r^2 * dr = 2*pi * spec * r^2 * dr
                E_contrib = 2.0_real64 * pi * spec(l, j) * r_j**2 * dr
                if (E_contrib == 0.0_real64) cycle

                K_tot = sqrt(real(l * (l + 1), real64)) / r_j

                ! Guard against log(0) for the DC / zero-wavenumber mode (l=0)
                if (K_tot <= 0.0_real64) then
                    E_K(0) = E_K(0) + E_contrib
                    cycle
                end if

                real_idx = (log(K_tot) - log_K_min) * inv_delta_log + 0.5_real64
                idx      = int(real_idx)

                if (idx < 1) then
                    E_K(1) = E_K(1) + E_contrib
                else if (idx >= n_bins) then
                    E_K(n_bins) = E_K(n_bins) + E_contrib
                else
                    w_hi         = real_idx - real(idx, real64)
                    w_lo         = 1.0_real64 - w_hi
                    E_K(idx)     = E_K(idx)     + E_contrib * w_lo
                    E_K(idx + 1) = E_K(idx + 1) + E_contrib * w_hi
                end if
            end do
        end do

        ! ---- Normalize by bin width (dK) to get spectral density E(K) ----
        w_lo = exp(0.5_real64 * delta_log) - exp(-0.5_real64 * delta_log)
        do j = 1, n_bins
            E_K(j) = E_K(j) / (K_bins(j) * w_lo)
        end do
    end subroutine bin_horizontal_spectrum

    ! =========================================================================
    ! VOLUME SPECTRUM ROUTINES
    ! =========================================================================

    ! 4. Computes spectral energy density per degree and shell for ONE component.
    !    energy(l, j) = sum_{m=0}^{l} [ C(l,m,j)^2 + S(l,m,j)^2 ]
    subroutine compute_spectral_energy(coeffs, lmax, n_shells, energy)
        integer, intent(in) :: lmax, n_shells
        real(kind=real64), intent(in)  :: coeffs(2, 0:lmax, 0:lmax, n_shells)
        real(kind=real64), intent(out) :: energy(0:lmax, n_shells)

        integer :: l, m, j

        energy = 0.0_real64
        do j = 1, n_shells
            do l = 0, lmax
                do m = 0, l
                    energy(l, j) = energy(l, j) &
                        + coeffs(1, l, m, j)**2 + coeffs(2, l, m, j)**2
                end do
            end do
        end do
    end subroutine compute_spectral_energy

    ! 5. Applies a 1-D Gaussian low-pass filter along the radial (shell) direction.
    !    Gaussian sigma = ell / sqrt(12), truncated at 3*sigma.
    !    Boundaries: truncate kernel and renormalize.
    subroutine radial_gaussian_filter(coeffs_in, lmax, n_shells, ell, coeffs_out)
        integer, intent(in) :: lmax, n_shells
        real(kind=real64), intent(in) :: ell
        real(kind=real64), intent(in)  :: coeffs_in(2, 0:lmax, 0:lmax, n_shells)
        real(kind=real64), intent(out) :: coeffs_out(2, 0:lmax, 0:lmax, n_shells)

        real(kind=real64) :: sigma, inv_2sig2
        integer :: j_max           ! half-width of the kernel in grid points
        real(kind=real64), allocatable :: w(:)   ! kernel weights, 0:2*j_max
        real(kind=real64) :: norm
        integer :: i, j, j_lo, j_hi, k

        if (ell <= 1.0e-8_real64) then
            coeffs_out = coeffs_in
            return
        end if

        sigma = ell / sqrt(12.0_real64)
        j_max = min(3 * ceiling(sigma), n_shells - 1)
        inv_2sig2 = 1.0_real64 / (2.0_real64 * sigma * sigma)

        ! Pre-compute symmetric kernel weights  w(0) .. w(2*j_max)
        !   w(k) = exp( -(k - j_max)^2 / (2 sigma^2) )
        allocate(w(0:2*j_max))
        do j = 0, 2*j_max
            w(j) = exp(-real((j - j_max)**2, real64) * inv_2sig2)
        end do

        ! Convolve each coefficient independently along the shell (j) direction
        coeffs_out = 0.0_real64
        do i = 1, n_shells
            ! Clamp the kernel window to stay within [1, n_shells]
            j_lo = max(0, j_max - (i - 1))       ! first valid kernel index
            j_hi = min(2*j_max, j_max + (n_shells - i))  ! last valid kernel index

            ! Accumulate weighted sum and normalisation
            norm = 0.0_real64
            do j = j_lo, j_hi
                k = i + (j - j_max)             ! shell index being sampled
                coeffs_out(:, :, :, i) = coeffs_out(:, :, :, i) &
                    + w(j) * coeffs_in(:, :, :, k)
                norm = norm + w(j)
            end do
            coeffs_out(:, :, :, i) = coeffs_out(:, :, :, i) / norm
        end do

        deallocate(w)
    end subroutine radial_gaussian_filter

    ! 6. Aggregates components of the 3D energy map E_modes(0:lmax, ell, j, c)
    subroutine accumulate_energy(se_diff, lmax, n_shells, r_min, dr, k, comp, E_modes, parseval_lhs)
        integer,           intent(in)     :: lmax, n_shells, k, comp
        real(kind=real64), intent(in)     :: se_diff(0:lmax, n_shells)
        real(kind=real64), intent(in)     :: r_min, dr
        real(kind=real64), intent(inout)  :: E_modes(0:lmax, 0:n_shells, n_shells, 3)
        real(kind=real64), intent(inout)  :: parseval_lhs

        real(kind=real64) :: r_j, E_dens
        integer :: l, j

        do j = 1, n_shells
            r_j = r_min + real(j - 1, real64) * dr
            do l = 0, lmax
                E_dens = 2.0_real64 * pi * se_diff(l, j) * r_j**2 * dr
                E_modes(l, k, j, comp) = E_modes(l, k, j, comp) + E_dens
                parseval_lhs = parseval_lhs + E_dens
            end do
        end do
    end subroutine accumulate_energy

    ! 7. Computes the 3D energy map E_modes(0:lmax, ell, j, c) by sweeping Gaussian
    !    filter scales over all velocity components. E_modes(0:lmax, ell, j, c) is
    !    the total kinetic energy contribution from angular degree l, radial
    !    filter scale ell, shell j, and velocity component comp.
    subroutine compute_volume_spectrum(coeffs_all, lmax, n_shells, r_min, dr, E_modes)
        integer,           intent(in)  :: lmax, n_shells
        real(kind=real64), intent(in)  :: coeffs_all(2, 0:lmax, 0:lmax, n_shells, 3)
        real(kind=real64), intent(in)  :: r_min, dr
        ! E_modes(0:lmax, ell, j, comp): energy at degree l, filter scale ell, shell j.
        ! Dimensions: (0:lmax) x (0:n_shells) x n_shells
        real(kind=real64), intent(out) :: E_modes(0:lmax, 0:n_shells, n_shells, 3)

        real(kind=real64), allocatable :: coeffs_filt(:,:,:,:), coeffs_avg(:,:,:,:)
        ! Filter Scale Energy (fse) and Higher Scale Energy (hse)
        real(kind=real64), allocatable :: fse(:,:), hse(:,:)

        real(kind=real64) :: r_j, parseval_lhs, parseval_rhs, ell
        integer :: comp, l, j, k

        allocate(coeffs_avg(2, 0:lmax, 0:lmax, n_shells))
        allocate(coeffs_filt(2, 0:lmax, 0:lmax, n_shells))
        allocate(hse(0:lmax, n_shells), fse(0:lmax, n_shells))

        E_modes      = 0.0_real64
        parseval_lhs = 0.0_real64
        parseval_rhs = 0.0_real64

        ! ==== Loop over velocity components ====
        do comp = 1, 3

            ! k = 0 means infinite filtering scale, i.e. simple averaging
            coeffs_avg = 0.0_real64
            do j = 1, n_shells
                coeffs_avg(:,:,:,1) = coeffs_avg(:,:,:,1) + coeffs_all(:,:,:,j,comp)
            end do
            coeffs_avg(:,:,:,1) = coeffs_avg(:,:,:,1) / n_shells
            do j = 2, n_shells
                coeffs_avg(:,:,:,j) = coeffs_avg(:,:,:,1)
            end do
            call compute_spectral_energy(coeffs_avg(:,:,:,:), lmax, n_shells, fse)

            ! Accumulate energy contribution into the 3D map (k = 0)
            call accumulate_energy(fse, lmax, n_shells, r_min, dr, 0, comp, E_modes, parseval_lhs)

            ! Set hse: first iteration (k = 1) differentiates against k = 0 state
            hse = fse

            ! ---- Sweep over radial filter scales ----
            ! Field filtered at scale ell has cumulative energy of all scales < ell.
            ! Pair-wise differences represents energy between the 2 filtering scales.
            do k = 1, n_shells-1
                ! ell = N/2 / sqrt(k(k+1))    : geometric average between N/2k and N/2(k+1)
                ell = 0.5_real64 * real(n_shells, real64) / sqrt(real(k*(k+1), real64))

                call radial_gaussian_filter(coeffs_all(:,:,:,:,comp), lmax, n_shells, &
                                            ell, coeffs_filt)
                call compute_spectral_energy(coeffs_filt, lmax, n_shells, fse)
                call accumulate_energy(fse-hse, lmax, n_shells, r_min, dr, k, comp, &
                                       E_modes, parseval_lhs)
                hse = fse

                if (mod(k, 25) == 0) write(*,'(A,I4,A,I4,A,I2)') &
                    "      k = ", k, " * 2*pi/(", n_shells, "dr),   component ", comp
            end do   ! ek

            ! Unfiltered spectral energy (k_max => ell = 0)
            call compute_spectral_energy(coeffs_all(:,:,:,:,comp), lmax, n_shells, fse)
            call accumulate_energy(fse-hse, lmax, n_shells, r_min, dr, n_shells, comp, E_modes, parseval_lhs)

            ! Parseval RHS: total energy from original field (last fse = unfiltered)
            do j = 1, n_shells
                r_j = r_min + real(j - 1, real64) * dr
                do l = 0, lmax
                    parseval_rhs = parseval_rhs + 2.0_real64 * pi * fse(l, j) * r_j**2 * dr
                end do
            end do

        end do   ! comp

        ! ---- Parseval diagnostic ----
        write(*,'(A)')          "   --- Parseval energy-conservation check ---"
        write(*,'(A, ES14.7)')  "       Binned spectrum sum  : ", parseval_lhs
        write(*,'(A, ES14.7)')  "       Direct original integral  : ", parseval_rhs
        if (parseval_rhs > 0.0_real64) then
            write(*,'(A, F8.4, A)') "       Relative difference: ", &
                abs(parseval_lhs - parseval_rhs) / parseval_rhs * 100.0_real64, " %"
        end if

        deallocate(coeffs_filt, hse, fse, coeffs_avg)
    end subroutine compute_volume_spectrum

    ! 8. Bins the 3D energy map E_modes(0:lmax, ell, j) into n_bins logarithmically-
    !    spaced bins.  Bin edges are determined from the theoretical min/max K values.
    !    Uses CIC interpolation between the two nearest bin centres.
    !    mode_flag: 0=all (default), 1=2D flow only, 2=Rest only.
    subroutine bin_spectrum(E_modes, lmax, n_shells, r_min, r_max, n_bins, E_K, K_bins, mode_flag)
        integer,           intent(in)  :: lmax, n_shells, n_bins
        real(kind=real64), intent(in)  :: E_modes(0:lmax, 0:n_shells, n_shells, 3)
        real(kind=real64), intent(in)  :: r_min, r_max
        real(kind=real64), intent(out) :: E_K(0:n_bins), K_bins(0:n_bins)
        integer, intent(in), optional  :: mode_flag

        real(kind=real64) :: dr, log_K_min, log_K_max, inv_delta_log, delta_log
        real(kind=real64) :: r_j, k_ang, k_rad, K_tot, E_contrib, real_idx, w_lo, w_hi
        integer :: l, j, k, comp, idx, flag

        ! ! DEBUG ONLY: Local arrays to count elements per bin (remove later)
        ! real(kind=real64) :: count_K(0:n_bins)
        ! ! DEBUG ONLY: 2D table of K_tot for each (radius j, filter k, degree l) (remove later)
        ! real(kind=real64), allocatable :: K_tot_table(:,:)
        ! integer :: debug_unit
        ! character(len=64) :: debug_fname

        flag = 0
        if (present(mode_flag)) flag = mode_flag

        ! count_K = 0.0_real64 ! DEBUG ONLY

        dr = (r_max - r_min) / real(n_shells - 1, real64)

        ! ---- Determine bin range from actual K values in the data ----
        ! K_min must be strictly positive since log(0) or log(negative) is undefined.
        ! We set K_min to the smallest non-zero physical scale (fundamental mode).


        log_K_min     = log(min(sqrt(2.0_real64)/r_max, 2*pi/(r_max-r_min)))
        log_K_max     = 0.5_real64 * log(real(lmax*(lmax+1), real64)/r_min**2 &
                      + (2*sqrt(real(n_shells*(n_shells+1), real64))*2*pi/(r_max-r_min))**2)
        delta_log     = (log_K_max - log_K_min) / real(n_bins, real64)
        inv_delta_log = real(n_bins, real64) / (log_K_max - log_K_min)

        ! ---- Build logarithmic bin centres directly ----
        K_bins(0) = 0.0_real64
        do j = 1, n_bins
            K_bins(j) = exp(log_K_min + (real(j, real64) - 0.5_real64) * delta_log)
        end do

        E_K = 0.0_real64

        ! ! DEBUG ONLY: Allocate K_tot table (rows = (j,l) pairs, cols = k values) (remove later)
        ! allocate(K_tot_table(n_shells * (lmax + 1), 0:n_shells))
        ! K_tot_table = 0.0_real64

        ! ---- CIC accumulation ----
        do comp = 1, 3
            ! mode_flag = 1: 2D flow (only horizontal components vt, vp)
            if (flag == 1 .and. comp == 2) cycle
            do j = 1, n_shells
                r_j = r_min + real(j - 1, real64) * dr
                do k = 0, n_shells
                    ! mode_flag = 1: 2D flow (only k=0 in horizontal components)
                    if (flag == 1 .and. k /= 0) cycle
                    ! mode_flag = 2: 3D flow (everything except 2D flow at k=0)
                    if (flag == 2 .and. (k == 0 .and. comp /= 2)) cycle

                    do l = 0, lmax
                        E_contrib = E_modes(l, k, j, comp)
                        k_ang = sqrt(real(l*(l+1), real64)) / r_j
                        k_rad = 2 * sqrt(real(k*(k+1), real64)) * 2*pi/(r_max-r_min)
                        K_tot = sqrt(k_ang**2 + k_rad**2)

                        ! ! DEBUG ONLY: Capture K_tot for comp==1 and flag==0 (remove later)
                        ! if (flag == 0 .and. comp == 1) &
                        !     K_tot_table((j-1)*(lmax+1) + l + 1, k) = K_tot

                        ! Guard against log(0) for the DC / zero-wavenumber mode (k=0, l=0)
                        if (K_tot <= 0.0_real64) then
                            E_K(0) = E_K(0) + E_contrib
                            ! count_K(0) = count_K(0) + 1.0_real64 ! DEBUG ONLY
                            cycle
                        end if

                        ! Direct O(1) index relative to bin centres (offset by +0.5)
                        real_idx = (log(K_tot) - log_K_min) * inv_delta_log + 0.5_real64
                        idx      = int(real_idx)

                        ! CIC: interpolate in logarithmic scale between the two nearest bin centres
                        if (idx < 1) then
                            E_K(1) = E_K(1) + E_contrib
                            ! count_K(1) = count_K(1) + 1.0_real64 ! DEBUG ONLY
                        else if (idx >= n_bins) then
                            E_K(n_bins) = E_K(n_bins) + E_contrib
                            ! count_K(n_bins) = count_K(n_bins) + 1.0_real64 ! DEBUG ONLY
                        else
                            w_hi         = real_idx - real(idx, real64)
                            w_lo         = 1.0_real64 - w_hi
                            E_K(idx)     = E_K(idx)     + E_contrib * w_lo
                            E_K(idx + 1) = E_K(idx + 1) + E_contrib * w_hi
                            ! count_K(idx)     = count_K(idx)     + w_lo ! DEBUG ONLY
                            ! count_K(idx + 1) = count_K(idx + 1) + w_hi ! DEBUG ONLY
                        end if
                    end do
                end do
            end do
        end do   ! comp

        ! ---- Normalize by bin width (dK) to get spectral density E(K) ----
        w_lo = exp(0.5_real64 * delta_log) - exp(-0.5_real64 * delta_log)
        do j = 1, n_bins
            E_K(j) = E_K(j) / (K_bins(j) * w_lo)
            ! count_K(j) = count_K(j) / (K_bins(j) * w_lo) ! DEBUG ONLY
        end do

        ! ! DEBUG ONLY: Write mode counts to file directly (remove later)
        ! debug_unit = 98
        ! if (flag == 0) then
        !     debug_fname = 'binning_count.dat'
        ! else if (flag == 1) then
        !     debug_fname = 'binning_count_2d.dat'
        ! else
        !     debug_fname = 'binning_count_3d.dat'
        ! end if
        ! open(unit=debug_unit, file=trim(debug_fname), status='replace')
        ! write(debug_unit, '(A)') "# k | normalized_count"
        ! do j = 0, n_bins
        !     write(debug_unit, '(ES15.7, ES15.7)') K_bins(j), count_K(j)
        ! end do
        ! close(debug_unit)

        ! ! DEBUG ONLY: Write K_tot table to file and deallocate (remove later)
        ! if (flag == 0) then
        !     debug_unit = 99
        !     open(unit=debug_unit, file='ktot_table.dat', status='replace')
        !     write(debug_unit, '(A)') '# K_tot table: rows=(j,l) pairs j=1..N_shells l=0..lmax, cols=k values k=0..N_shells'
        !     do j = 1, n_shells * (lmax + 1)
        !         write(debug_unit, '(*(ES12.4, 1X))') K_tot_table(j, :)
        !     end do
        !     close(debug_unit)
        ! end if
        ! deallocate(K_tot_table)
    end subroutine bin_spectrum

    ! =========================================================================
    ! OGIVE (CUMULATIVE-SPLINE) SPECTRUM ROUTINES
    ! =========================================================================

    ! 9. Recursive quicksort on key(:), permuting vals(:) concurrently.
    recursive subroutine quicksort_paired(key, vals, lo_in, hi_in)
        real(kind=real64), intent(inout) :: key(:), vals(:)
        integer, intent(in) :: lo_in, hi_in

        real(kind=real64) :: pivot, tmp
        integer :: lo, hi, mid

        lo = lo_in
        hi = hi_in
        if (lo >= hi) return

        ! Median-of-three pivot selection
        mid = (lo + hi) / 2
        if (key(lo) > key(mid)) then
            tmp = key(lo); key(lo) = key(mid); key(mid) = tmp
            tmp = vals(lo); vals(lo) = vals(mid); vals(mid) = tmp
        end if
        if (key(lo) > key(hi)) then
            tmp = key(lo); key(lo) = key(hi); key(hi) = tmp
            tmp = vals(lo); vals(lo) = vals(hi); vals(hi) = tmp
        end if
        if (key(mid) > key(hi)) then
            tmp = key(mid); key(mid) = key(hi); key(hi) = tmp
            tmp = vals(mid); vals(mid) = vals(hi); vals(hi) = tmp
        end if
        pivot = key(mid)

        ! Hoare partition
        lo = lo_in
        hi = hi_in
        do
            do while (key(lo) < pivot)
                lo = lo + 1
            end do
            do while (key(hi) > pivot)
                hi = hi - 1
            end do
            if (lo >= hi) exit
            tmp = key(lo); key(lo) = key(hi); key(hi) = tmp
            tmp = vals(lo); vals(lo) = vals(hi); vals(hi) = tmp
            lo = lo + 1
            hi = hi - 1
        end do

        call quicksort_paired(key, vals, lo_in, hi)
        call quicksort_paired(key, vals, hi + 1, hi_in)
    end subroutine quicksort_paired

    ! 9b. Merge duplicate wavenumbers in place after sorting.
    !     Returns the number of unique modes in n_out.
    !     Optionally zeros out the remaining elements for safety/debugging.
    subroutine merge_duplicates(key, vals, n_in, n_out)
        real(kind=real64), intent(inout) :: key(:), vals(:)
        integer, intent(in)  :: n_in
        integer, intent(out) :: n_out

        integer :: i, j

        if (n_in <= 1) then
            n_out = n_in
            return
        end if

        j = 1  ! write pointer
        do i = 2, n_in
            if (key(i) == key(j)) then
                ! Duplicate found: sum energy
                vals(j) = vals(j) + vals(i)
            else
                ! New unique key: advance and copy
                j = j + 1
                key(j) = key(i)
                vals(j) = vals(i)
            end if
        end do
        n_out = j

        ! Zero out the remaining unused elements
        if (n_out < n_in) then
            key(n_out+1:n_in) = 0.0_real64
            vals(n_out+1:n_in) = 0.0_real64
        end if
    end subroutine merge_duplicates

    ! 10. Ogive binning for the flattened energy E_flat.
    !     Sorts, cumulates, fits a smoothing spline, and analytically
    !     differentiates onto a clean log-spaced target grid.
    subroutine ogive_binning(E_flat, K_flat, n_modes, K_min, K_max, n_bins, E_K, K_bins)
        integer,           intent(in)  :: n_modes, n_bins
        real(kind=real64), intent(inout) :: E_flat(n_modes), K_flat(n_modes)
        real(kind=real64), intent(in)  :: K_min, K_max
        real(kind=real64), intent(out) :: E_K(0:n_bins), K_bins(0:n_bins)

        real(kind=real64), allocatable :: Cumul(:)
        real(kind=real64) :: dK, s_factor
        integer :: idx, n_unique

        type(fitpack_curve) :: spl
        integer(FP_FLAG) :: ierr


        ! ---- Step 1: Sort by wavenumber (quicksort) ----
        call quicksort_paired(K_flat, E_flat, 1, n_modes)

        ! ---- Step 1b: Merge duplicate wavenumbers ----
        call merge_duplicates(K_flat, E_flat, n_modes, n_unique)

        ! ---- Step 2: Cumulative sum ----
        allocate(Cumul(n_unique))
        Cumul(1) = E_flat(1)
        do idx = 2, n_unique
            Cumul(idx) = Cumul(idx - 1) + E_flat(idx)
        end do

        write(*,'(A, ES12.5)') "   Ogive: total cumulative energy = ", Cumul(n_unique)

        s_factor = s_factor_mult * n_unique
        ierr = spl%new_fit(K_flat(1:n_unique), Cumul, smoothing=s_factor)
        if (.not. FITPACK_SUCCESS(ierr)) then
            write(error_unit, '(A, I0)') "   Ogive: spline fit error, code = ", ierr
            write(error_unit, '(A)')     "   Ogive: " // FITPACK_MESSAGE(ierr)
            E_K = 0.0_real64
            K_bins = 0.0_real64
            deallocate(Cumul)
            return
        end if

        ! ---- Step 4: Generate clean log-spaced target grid ----
        ! Use the theoretical K_max so the target grid is identical across all calls
        dK = exp((log(K_max) - log(K_min)) / real(n_bins - 1, real64))
        K_bins(0) = 0.0_real64
        K_bins(1) = K_min
        do idx = 2, n_bins
            K_bins(idx) = K_bins(idx - 1) * dK
        end do

        ! ---- Step 5: Analytical differentiation dI/dK ----
        E_K = spl%dfdx(K_bins, 1, ierr)
        if (.not. FITPACK_SUCCESS(ierr)) then
            write(error_unit, '(A, I0)') "   Ogive: derivative evaluation error, code = ", ierr
            E_K(0:n_bins) = 0.0_real64
            deallocate(Cumul)
            return
        end if

        ! Zero out any extrapolation beyond the maximum actual data wavenumber
        where (K_bins > K_flat(n_unique)) E_K = 0.0_real64

        ! Clip tiny negative values that may arise from spline oscillation
        where (E_K < 0.0_real64) E_K = 0.0_real64

        deallocate(Cumul)
        call spl%destroy()
    end subroutine ogive_binning

    ! 11. Ogive spectrum for the 3D volume energy map E_modes.
    !     Gathers discrete modes, and returns a log-spaced grid.
    !     mode_flag: 0=all (default), 1=2D modes of 2D flow only,
    !      2=3D, all the rest.
    subroutine ogive_spectrum(E_modes, lmax, n_shells, r_min, r_max, &
                              n_bins, E_K, K_bins, mode_flag)
        integer,           intent(in)  :: lmax, n_shells, n_bins
        real(kind=real64), intent(in)  :: E_modes(0:lmax, 0:n_shells, n_shells, 3)
        real(kind=real64), intent(in)  :: r_min, r_max
        real(kind=real64), intent(out) :: E_K(0:n_bins), K_bins(0:n_bins)
        integer, intent(in), optional  :: mode_flag

        real(kind=real64) :: E_modes_cs(0:lmax, 0:n_shells, n_shells)
        real(kind=real64) :: dr, K_min, K_max, vol_fac
        real(kind=real64) :: r_j, k_ang, k_rad, K_tot, E_contrib
        integer :: l, j, k, comp, idx, n_modes, flag

        real(kind=real64), allocatable :: K_flat(:), E_flat(:)

        flag = 0
        if (present(mode_flag)) flag = mode_flag

        dr = (r_max - r_min) / real(n_shells - 1, real64)

        ! ---- Sum spacial components ----
        E_modes_cs = 0.0_real64
        do comp = 1, 3
            ! mode_flag = 1: 2D flow (only horizontal components vt, vp)
            if (flag == 1 .and. comp == 2) cycle
            E_modes_cs = E_modes_cs + E_modes(:,:,:, comp)
        end do

        ! ---- Count non-zero modes to pre-allocate ----
        n_modes = 0
        do j = 1, n_shells
            do k = 0, n_shells
                ! mode_flag = 1: 2D flow (only k=0 in horizontal components)
                if (flag == 1 .and. k /= 0) cycle
                do l = 0, lmax
                    ! mode_flag = 2: 3D flow (everything except 2D flow at k=0)
                    if (flag == 2 .and. k == 0) then
                        E_contrib = E_modes(l, k, j, 2)
                    else
                        E_contrib = E_modes_cs(l, k, j)
                    end if
                    if (E_contrib == 0.0_real64) cycle
                    n_modes = n_modes + 1
                end do
            end do
        end do

        write(*,'(A, I0, A)') "   Ogive: gathered ", n_modes, " non-zero modes"

        if (n_modes == 0) then
            E_K = 0.0_real64
            K_bins = 0.0_real64
            return
        end if

        ! ---- Gather modes into flat arrays ----
        allocate(K_flat(n_modes), E_flat(n_modes))
        vol_fac = 4 * pi/(r_max - r_min)
        idx = 0
        do j = 1, n_shells
            do k = 0, n_shells
                ! mode_flag = 1: 2D flow (only k=0 in horizontal components)
                if (flag == 1 .and. k /= 0) cycle
                do l = 0, lmax
                    ! mode_flag = 2: 3D flow (everything except 2D flow at k=0)
                    if (flag == 2 .and. k == 0) then
                        E_contrib = E_modes(l, k, j, 2)
                    else
                        E_contrib = E_modes_cs(l, k, j)
                    end if
                    if (E_contrib == 0.0_real64) cycle
                    r_j = r_min + real(j - 1, real64) * dr
                    k_ang = sqrt(real(l*(l+1), real64)) / r_j
                    k_rad = sqrt(real(k*(k+1), real64)) * vol_fac
                    K_tot = sqrt(k_ang**2 + k_rad**2)
                    idx = idx + 1
                    K_flat(idx) = K_tot
                    E_flat(idx) = E_contrib
                end do
            end do
        end do

        K_min = min(sqrt(2.0_real64)/r_max, 2*pi/(r_max-r_min))
        K_max = real(lmax*(lmax+1), real64)/r_min**2 &
              + (sqrt(real(n_shells*(n_shells+1), real64))*vol_fac)**2

        call ogive_binning(E_flat, K_flat, n_modes, K_min, K_max, n_bins, E_K, K_bins)

        deallocate(K_flat, E_flat)

    end subroutine ogive_spectrum

    ! 12. Ogive spectrum for the 2D horizontal spectrum spec.
    !     Gathers discrete modes, and returns a log-spaced grid.
    subroutine ogive_horizontal_spectrum(spec, lmax, n_shells, r_min, r_max, &
                                         n_bins, E_K, K_bins)
        integer,           intent(in)  :: lmax, n_shells, n_bins
        real(kind=real64), intent(in)  :: spec(0:lmax, n_shells)
        real(kind=real64), intent(in)  :: r_min, r_max
        real(kind=real64), intent(out) :: E_K(0:n_bins), K_bins(0:n_bins)

        real(kind=real64) :: dr, K_min, K_max
        real(kind=real64) :: r_j, K, E_contrib
        integer :: l, j, idx, n_modes

        real(kind=real64), allocatable :: K_flat(:), E_flat(:)

        dr = (r_max - r_min) / real(n_shells - 1, real64)

        ! ---- Count non-zero modes to pre-allocate ----
        n_modes = 0
        do j = 1, n_shells
            r_j = r_min + real(j - 1, real64) * dr
            if (r_j == 0.0_real64) cycle
            do l = 0, lmax
                if (spec(l, j) == 0.0_real64) cycle
                n_modes = n_modes + 1
            end do
        end do

        write(*,'(A, I0, A)') "   Horizontal Ogive: gathered ", n_modes, " non-zero modes"

        if (n_modes == 0) then
            E_K = 0.0_real64
            K_bins = 0.0_real64
            return
        end if

        ! ---- Gather modes into flat arrays ----
        allocate(K_flat(n_modes), E_flat(n_modes))
        idx = 0
        do j = 1, n_shells
            r_j = r_min + real(j - 1, real64) * dr
            do l = 0, lmax
                E_contrib = 2.0_real64 * pi * spec(l, j) * r_j**2 * dr
                if (E_contrib == 0.0_real64) cycle
                K = sqrt(real(l*(l+1), real64)) / r_j
                idx = idx + 1
                K_flat(idx) = K
                E_flat(idx) = E_contrib
            end do
        end do

        K_min = min(sqrt(2.0_real64)/r_max, 2*pi/(r_max-r_min))
        K_max = real(lmax*(lmax+1), real64)/r_min**2 &
              + (2*sqrt(real(n_shells*(n_shells+1), real64))*2*pi/(r_max-r_min))**2

        call ogive_binning(E_flat, K_flat, n_modes, K_min, K_max, n_bins, E_K, K_bins)

        deallocate(K_flat, E_flat)

    end subroutine ogive_horizontal_spectrum

end module sht_processing_mod


program sht_3d_utility
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    use data_io_mod
    use sht_processing_mod
    implicit none

    integer :: n1, n2, n3, lmax, n_lat, n_lon
    integer :: arg_count, i, j, k, io_status
    real(kind=real64) :: sim_time, pi, theta, phi, ct, st, cp, sp

    character(len=256) :: input_path, coord_type, out_file, r_arg, arg
    integer :: arg_idx
    logical :: r_found
    integer, parameter :: out_unit = 11, vol_unit = 12

    real(kind=real64), allocatable :: v1(:,:,:), v2(:,:,:), v3(:,:,:)
    ! Velocity in spherical coordinates (vt=azimuthal/zonal, vr=radial, vp=polar/meridional)
    real(kind=real64), allocatable :: vt(:,:,:),vr(:,:,:),  vp(:,:,:)

    real(kind=real64), allocatable :: grid_t(:,:), grid_r(:,:), grid_p(:,:)
    real(kind=real64), allocatable :: coeffs_t(:,:,:), coeffs_r(:,:,:), coeffs_p(:,:,:)
    real(kind=real64), allocatable :: spec_tot(:,:), spec_m0(:,:), spec_rest(:,:)

    ! Volume spectrum arrays
    real(kind=real64), allocatable :: coeffs_all(:,:,:,:,:)
    real(kind=real64), allocatable :: E_modes(:,:,:,:)
    real(kind=real64), allocatable :: E_K(:), E2d_K(:), E3d_K(:), K_bins(:)
    real(kind=real64), allocatable :: E_tot_H(:), E_zonal_H(:), E_rest_H(:)
    real(kind=real64) :: rint, rext, dr
    integer :: n_bins

    real(kind=real64) :: vx, vy, vz
    real(kind=real64) :: E_phys, d_theta, d_phi, r_j
    real(kind=real64) :: E_binned_sum!, log_K_min, log_K_max, delta_log, w_hi

    ! 1. Parse Command Line
    arg_count = command_argument_count()
    if (arg_count < 2) then
        write(error_unit, '("Usage: ./sht_3d_utility <file_path> -r=<r_in> [-s=<s_factor>] [-p|-c]")')
        write(error_unit, '(" -r=<r_in>      : inner radius (required); r_ext = r_in + 1")')
        write(error_unit, '(" -s=<s_factor>  : smoothing factor multiplier (default=0.005)")')
        write(error_unit, '(" -p (default)   = velocity field input is in spherical coordinates")')
        write(error_unit, '(" -c             = velocity field input is in Cartesian coordinates")')
        write(error_unit, '("Example: ./sht_3d_utility dati/spbou001234.q -r=1.0 -s=0.002 -c")')
        stop 1
    end if

    call get_command_argument(1, input_path)

    ! Default values
    coord_type = '-p'
    r_found = .false.

    do arg_idx = 2, arg_count
        call get_command_argument(arg_idx, arg)
        if (arg(1:3) == '-r=') then
            read(arg(4:), *, iostat=io_status) rint
            if (io_status /= 0) then
                write(error_unit, '("Error: could not parse r_in from: ", A)') trim(arg)
                stop 1
            end if
            r_found = .true.
        else if (arg(1:3) == '-s=') then
            read(arg(4:), *, iostat=io_status) s_factor_mult
            if (io_status /= 0) then
                write(error_unit, '("Error: could not parse s_factor from: ", A)') trim(arg)
                stop 1
            end if
        else if (trim(arg) == '-p' .or. trim(arg) == '-c') then
            coord_type = arg
        else
            write(error_unit, '("Warning: unknown argument skipped: ", A)') trim(arg)
        end if
    end do

    if (.not. r_found) then
        write(error_unit, '("Error: inner radius -r=<r_in> is required.")')
        stop 1
    end if

    rext = rint + 1.0_real64
    write(*, '(A, F10.4, A, F10.4)') "   Inner radius: ", rint, "   Outer radius: ", rext
    write(*, '(A, F10.6)') "   Ogive smoothing factor multiplier: ", s_factor_mult

    ! 2. Read Data
    call read_simulation_data(input_path, n1, n2, n3, sim_time, v1, v2, v3)
    write(*, '(A, 3(I0, A), F10.4)') "   Dimensions: ", n1, "x", n2, "x", n3, " at time ", sim_time

    ! Compute radial spacing (uniform grid, shell thickness = 1)
    dr = (rext - rint) / real(n2 - 1, real64)

    ! 3. Handle Velocity Field Appropriately (Polar vs Cartesian)
    allocate(vt(n1, n2, n3), vr(n1, n2, n3), vp(n1, n2, n3))

    if (trim(coord_type) == '-p') then
        write(*,*) "-> Using spherical velocity components directly (vt, vr, vp)."
        vt = v1;   vr = v2;   vp = v3
    else if (trim(coord_type) == '-c') then
        write(*,*) "-> Converting Cartesian velocity components to spherical..."
        pi = 2.0_real64 * asin(1.0_real64)

        do k = 1, n3
            ! Colatitude (polar angle phi) varies along k
            phi = real(k-1, real64) * pi / real(n3-1, real64)
            cp = cos(phi)
            sp = sin(phi)
            do j = 1, n2
                do i = 1, n1
                    ! Longitude (azimuthal angle theta) varies along i
                    theta = real(i-1, real64) * 2.0_real64 * pi / real(n1-1, real64)
                    ct = cos(theta)
                    st = sin(theta)

                    vx = v1(i, j, k)
                    vy = v2(i, j, k)
                    vz = v3(i, j, k)

                    vt(i, j, k) = -vx*st    + vy*ct
                    vr(i, j, k) =  vx*sp*ct + vy*sp*st + vz*cp
                    vp(i, j, k) =  vx*cp*ct + vy*cp*st - vz*sp
                end do
            end do
        end do
    else
        write(error_unit, '("Error: Unknown coordinate type option. Must be ''-p'' or ''-c''.")')
        stop 1
    end if

    ! Calculate total physical energy from the spatial grid
    E_phys = 0.0_real64
    d_phi = pi / real(n3-1, real64)
    d_theta = 2.0_real64 * pi / real(n1-1, real64)
    do k = 1, n3
        phi = real(k-1, real64) * pi / real(n3-1, real64)
        sp = sin(phi)
        do j = 1, n2
            r_j = rint + real(j-1, real64) * dr
            ! Only sum i=1 to n1-1 since i=n1 is a periodic copy of i=1
            do i = 1, n1-1
                E_phys = E_phys + 0.5_real64 * (vt(i,j,k)**2 + vr(i,j,k)**2 + vp(i,j,k)**2) * &
                         r_j**2 * dr * d_theta * d_phi * sp
            end do
        end do
    end do
    write(*,'(A, ES14.7)') "   Total Physical Energy (Grid sum) : ", E_phys

    ! We don't need the original arrays anymore
    deallocate(v1, v2, v3)

    ! 4. Setup SHT Parameters
    n_lat = n3 - mod(n3, 2)
    n_lon = n1 - 1
    lmax = n_lat / 2 - 1

    allocate(grid_t(n_lat, n_lon), grid_r(n_lat, n_lon), grid_p(n_lat, n_lon))
    allocate(coeffs_t(2, 0:lmax, 0:lmax), coeffs_r(2, 0:lmax, 0:lmax), coeffs_p(2, 0:lmax, 0:lmax))

    allocate(spec_tot(0:lmax, n2), spec_m0(0:lmax, n2), spec_rest(0:lmax, n2))

    ! Allocate storage for all SH coefficients (all shells, all 3 components)
    allocate(coeffs_all(2, 0:lmax, 0:lmax, n2, 3))

    spec_tot = 0.0_real64
    spec_m0 = 0.0_real64
    spec_rest = 0.0_real64

    write(*,*) "-> Performing SHT on Velocity Field..."
    do j = 1, n2
        ! Extract 2D spherical slices, transpose to (lat, lon) format for SHTOOLS
        grid_t(:, :) = transpose(vt(1:n_lon, j, 1:n_lat))
        grid_r(:, :) = transpose(vr(1:n_lon, j, 1:n_lat))
        grid_p(:, :) = transpose(vp(1:n_lon, j, 1:n_lat))

        ! Call Subroutine 1: Get Coefficients
        call get_sht_coeffs(grid_t, n_lat, lmax, coeffs_t)
        call get_sht_coeffs(grid_r, n_lat, lmax, coeffs_r)
        call get_sht_coeffs(grid_p, n_lat, lmax, coeffs_p)

        ! Store coefficients for volume spectrum
        coeffs_all(:, :, :, j, 1) = coeffs_t
        coeffs_all(:, :, :, j, 2) = coeffs_r
        coeffs_all(:, :, :, j, 3) = coeffs_p

        ! Call Subroutine 2: Get Energy Spectra
        call upd_energy_spectrum(coeffs_t, lmax, spec_tot(:, j), &
                                 spec_m0(:, j), spec_rest(:, j))
        call upd_energy_spectrum(coeffs_r, lmax, spec_tot(:, j), &
                                 spec_m0(:, j), spec_rest(:, j))
        call upd_energy_spectrum(coeffs_p, lmax, spec_tot(:, j), &
                                 spec_m0(:, j), spec_rest(:, j))

        if (mod(j, 25) == 0) write(*,*) "      Progress: shell ", j, "/", n2
    end do

    deallocate(vt, vr, vp)
    deallocate(grid_t, grid_r, grid_p)
    deallocate(coeffs_t, coeffs_r, coeffs_p)

    ! 5. Bin the horizontal spectra
    ! Number of logarithmic bins (~12 per octave)
    n_bins = max(64, nint(12 * log(real(max(lmax, n2), real64)) / log(2.0_real64)))
    allocate(E_tot_H(0:n_bins), E_zonal_H(0:n_bins), E_rest_H(0:n_bins), K_bins(0:n_bins))

    write(*,*) "-> Computing Horizontal Spectra (ogive)..."
    call ogive_horizontal_spectrum(spec_tot,  lmax, n2, rint, rext, n_bins, E_tot_H,   K_bins)
    call ogive_horizontal_spectrum(spec_m0,   lmax, n2, rint, rext, n_bins, E_zonal_H, K_bins)
    call ogive_horizontal_spectrum(spec_rest, lmax, n2, rint, rext, n_bins, E_rest_H,  K_bins)

    ! 6. Output the horizontal spectra
    out_file = 'horizontal_spectrum.dat'
    write(*,*) "-> Saving Horizontal Spectra to: ", trim(out_file)
    open(unit=out_unit, file=trim(out_file), status='replace')
    write(out_unit, '(A)') "# k | E_tot_H(k) | E_zonal_H(k) | E_rest_H(k)"
    do i = 0, n_bins
        if (E_tot_H(i) > 0.0_real64) then
            write(out_unit, '(ES15.7, ES15.7, ES15.7, ES15.7)') &
                K_bins(i), E_tot_H(i), E_zonal_H(i), E_rest_H(i)
        end if
    end do
    close(out_unit)

    deallocate(E_tot_H, E_zonal_H, E_rest_H)
    deallocate(spec_tot, spec_m0, spec_rest)

    ! 7. Compute and Output the 3D Volume Spectrum E(K)
    write(*,*) "-> Computing 3D Volume Spectrum E(K)..."

    ! E_modes(0:lmax, ell, j): energy at degree l, filter scale ell, shell j, component
    allocate(E_modes(0:lmax, 0:n2, n2, 3))

    call compute_volume_spectrum(coeffs_all, lmax, n2, rint, dr, E_modes)
    deallocate(coeffs_all)

    allocate(E_K(0:n_bins), E2d_K(0:n_bins), E3d_K(0:n_bins))
    write(*,*) "-> Computing Volume Spectra..."
    call ogive_spectrum(E_modes, lmax, n2, rint, rext, n_bins, E_K, K_bins)
    call ogive_spectrum(E_modes, lmax, n2, rint, rext, n_bins, E2d_K, K_bins, 1)
    call ogive_spectrum(E_modes, lmax, n2, rint, rext, n_bins, E3d_K, K_bins, 2)

    ! ! Parseval check: integrate binned spectrum as sum(E_K * dK)
    ! delta_log = log(K_bins(2) / K_bins(1))
    ! w_hi = exp(0.5_real64 * delta_log) - exp(-0.5_real64 * delta_log)
    
    ! E_binned_sum = E_K(0) ! DC mode doesn't have a width in the log bins
    ! do i = 1, n_bins
    !     E_binned_sum = E_binned_sum + E_K(i) * K_bins(i) * w_hi
    ! end do

    ! Parseval check: integrate E(K) dK via trapezoidal rule on the ogive output
    
    E_binned_sum = 0.0_real64
    do i = 1, n_bins
        E_binned_sum = E_binned_sum &
            + 0.5_real64 * (E_K(i-1) + E_K(i)) * (K_bins(i) - K_bins(i-1))
    end do

    write(*,'(A)')          "   --- Final Physical Parseval Check ---"
    write(*,'(A, ES14.7)')  "       Physical Energy (Grid sum)   : ", E_phys
    write(*,'(A, ES14.7)')  "       Integral of Binned Spectrum  : ", E_binned_sum
    if (E_phys > 0.0_real64) then
        write(*,'(A, F8.4, A)') "       Relative difference: ", &
            abs(E_binned_sum - E_phys) / E_phys * 100.0_real64, " %"
    end if

    deallocate(E_modes)

    ! Write volume spectrum
    out_file = 'volume_spectrum.dat'
    write(*,*) "-> Saving Volume Spectrum to: ", trim(out_file)
    open(unit=vol_unit, file=trim(out_file), status='replace')
    write(vol_unit, '(A)') "# k | E(K) | E_2d(k) | E_3d(k)"
    do i = 0, n_bins
        if (E_K(i) > 0.0_real64) then
            write(vol_unit, '(ES15.7, ES15.7, ES15.7, ES15.7)') K_bins(i), E_K(i), E2d_K(i), E3d_K(i)
        end if
    end do
    close(vol_unit)

    deallocate(E_K, E2d_K, E3d_K, K_bins)

    write(*,*) "-> Done."

end program sht_3d_utility
