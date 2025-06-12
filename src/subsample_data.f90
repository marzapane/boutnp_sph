! =============================================================================
! Programma per leggere, sottocampionare e riscrivere dati da file binari
! generati da una simulazione Fortran 77.
!
! Scritto in Fortran 95/2008.
!
! Compilazione (con gfortran):
! gfortran -o src/subsampler subsample_data.f90 -Wall -Wextra
!
! Esecuzione:
! ./subsampler
! =============================================================================
program subsample_data
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  implicit none

  ! ---------------------------------------------------------------------------
  ! Variabili
  ! ---------------------------------------------------------------------------

  ! Dimensioni della griglia
  integer :: n1, n2, n3
  integer :: n1_check, n2_check, n3_check
  integer :: n1_new, n2_new, n3_new

  ! Variabili per i parametri fisici letti dal file .q
  real(kind=real64) :: epsil, re, sim_time

  ! Array allocabili per contenere i dati letti dai file originali
  real(kind=real64), allocatable, dimension(:,:,:) :: xxx, yyy, zzz
  real(kind=real64), allocatable, dimension(:,:,:) :: E, vx, vy, vz, pr

  ! Array allocabili per i dati sottocampionati da scrivere nei nuovi file
  real(kind=real64), allocatable, dimension(:,:,:) :: xxx_new, yyy_new, zzz_new
  real(kind=real64), allocatable, dimension(:,:,:) :: E_new, vx_new, vy_new, vz_new, pr_new

  ! Variabili per la gestione dei file e degli errori I/O
  character(len=256) :: grid_in_file, data_in_file
  character(len=256) :: grid_out_file, data_out_file
  integer            :: io_status, dummy_int
  integer, parameter :: in_unit = 10, out_unit = 11

  ! Stride per il sottocampionamento
  integer, parameter :: i_stride = 2, j_stride = 4, k_stride = 2

  ! ---------------------------------------------------------------------------
  ! Inizio del programma eseguibile
  ! ---------------------------------------------------------------------------

  ! 1. Ottiene i nomi dei file dall'utente
  write(*, '(a)', advance='no') "Inserisci nome del file della griglia di input (es. grid.xyz): "
  read(*, '(a)') grid_in_file
  write(*, '(a)', advance='no') "Inserisci nome del file dei dati di input (es. dati/boum000100.q): "
  read(*, '(a)') data_in_file
  write(*, *)
  write(*, '(a)', advance='no') "Inserisci nome del file della griglia di output (es. grid_sub.xyz): "
  read(*, '(a)') grid_out_file
  write(*, '(a)', advance='no') "Inserisci nome del file dei dati di output (es. data_sub.q): "
  read(*, '(a)') data_out_file
  write(*, *) "------------------------------------------------------------"

  ! ===========================================================================
  ! 2. Processa il file della griglia (grid.xyz)
  ! ===========================================================================
  write(*, '(a, a)') "-> Processando il file della griglia: ", trim(grid_in_file)
  open(unit=in_unit, file=trim(grid_in_file), form='unformatted', status='old', &
       action='read', iostat=io_status)
  if (io_status /= 0) stop "Errore: impossibile aprire il file di input della griglia."

  ! Lettura dei dati dal file grid.xyz secondo la struttura F77
  read(in_unit) dummy_int  ! Legge il primo record (intero '1')
  read(in_unit) n1, n2, n3  ! Legge il secondo record (dimensioni)

  write(*, '(a, 3(i0, a))') "   Dimensioni griglia originale: (", n1, ", ", n2, ", ", n3, ")"

  ! Allocazione degli array per i dati completi
  allocate(xxx(n1, n2, n3), yyy(n1, n2, n3), zzz(n1, n2, n3), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per la griglia completa."

  ! Lettura del quarto record contenente i 3 array delle coordinate
  read(in_unit) xxx, yyy, zzz

  close(in_unit)

  ! Calcolo delle nuove dimensioni dopo il sottocampionamento
  n1_new = size(xxx(1::i_stride, 1, 1))
  n2_new = size(xxx(1, 1::j_stride, 1))
  n3_new = size(xxx(1, 1, 1::k_stride))
  write(*, '(a, 3(i0, a))') "   Dimensioni griglia ridotta:   (", n1_new, ", ", n2_new, ", ", n3_new, ")"

  ! Allocazione degli array per i dati sottocampionati
  allocate(xxx_new(n1_new, n2_new, n3_new), &
           yyy_new(n1_new, n2_new, n3_new), &
           zzz_new(n1_new, n2_new, n3_new), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per la griglia ridotta."

  ! Esecuzione del sottocampionamento tramite "array slicing"
  xxx_new = xxx(1::i_stride, 1::j_stride, 1::k_stride)
  yyy_new = yyy(1::i_stride, 1::j_stride, 1::k_stride)
  zzz_new = zzz(1::i_stride, 1::j_stride, 1::k_stride)

  ! Deallocazione degli array originali non più necessari
  deallocate(xxx, yyy, zzz)

  ! Scrittura del nuovo file della griglia sottocampionato
  open(unit=out_unit, file=trim(grid_out_file), form='unformatted', status='replace', &
       action='write', iostat=io_status)
  if (io_status /= 0) stop "Errore: impossibile creare il file di output della griglia."
  write(out_unit) 1
  write(out_unit) n1_new, n2_new, n3_new
  write(out_unit) xxx_new, yyy_new, zzz_new
  close(out_unit)

  deallocate(xxx_new, yyy_new, zzz_new)

  write(*, '(a, a)') "   File griglia ridotta creato: ", trim(grid_out_file)
  write(*, *)

  ! ===========================================================================
  ! 3. Processa il file dei dati (.q)
  ! ===========================================================================
  write(*, '(a, a)') "-> Processando il file dei dati: ", trim(data_in_file)
  open(unit=in_unit, file=trim(data_in_file), form='unformatted', status='old', &
       action='read', iostat=io_status)
  if (io_status /= 0) stop "Errore: impossibile aprire il file di input dei dati."

  ! Lettura dell'header del file .q
  read(in_unit) dummy_int
  read(in_unit) n1_check, n2_check, n3_check

  ! Verifica che le dimensioni siano consistenti con il file della griglia
  if (n1_check /= n1 .or. n2_check /= n2 .or. n3_check /= n3) then
    write(error_unit, *) "Errore: Le dimensioni nel file .q non corrispondono a quelle del file grid.xyz."
    stop 2
  end if

  read(in_unit) epsil, epsil, re, sim_time

  ! Allocazione degli array per i dati completi
  allocate(E(n1, n2, n3), vx(n1, n2, n3), vy(n1, n2, n3), &
          vz(n1, n2, n3), pr(n1, n2, n3), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per i dati completi."

  ! Lettura dei 5 array dei dati
  read(in_unit) E, vx, vy, vz, pr

  close(in_unit)

  ! Allocazione degli array per i dati sottocampionati
  allocate(E_new(n1_new, n2_new, n3_new),  &
           vx_new(n1_new, n2_new, n3_new), &
           vy_new(n1_new, n2_new, n3_new), &
           vz_new(n1_new, n2_new, n3_new), &
           pr_new(n1_new, n2_new, n3_new), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per i dati ridotti."

  ! Sottocampionamento dei 5 array
  E_new  = E (1::i_stride, 1::j_stride, 1::k_stride)
  vx_new = vx(1::i_stride, 1::j_stride, 1::k_stride)
  vy_new = vy(1::i_stride, 1::j_stride, 1::k_stride)
  vz_new = vz(1::i_stride, 1::j_stride, 1::k_stride)
  pr_new = pr(1::i_stride, 1::j_stride, 1::k_stride)
  
  deallocate(E, vx, vy, vz, pr)

  ! Scrittura del nuovo file dei dati sottocampionato
  open(unit=out_unit, file=trim(data_out_file), form='unformatted', status='replace', &
       action='write', iostat=io_status)
  if (io_status /= 0) stop "Errore: impossibile creare il file di output dei dati."

  write(out_unit) 1
  write(out_unit) n1_new, n2_new, n3_new
  write(out_unit) epsil, epsil, re, sim_time ! Scrive i parametri originali
  write(out_unit) E_new, vx_new, vy_new, vz_new, pr_new

  close(out_unit)

  deallocate(E_new, vx_new, vy_new, vz_new, pr_new)

  write(*, '(a, a)') "   File dati ridotti creato: ", trim(data_out_file)

  write(*, *) "------------------------------------------------------------"
  write(*, *) "Operazione completata con successo."

end program subsample_data
