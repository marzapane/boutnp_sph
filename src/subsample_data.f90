! =============================================================================
! Programma per leggere, sottocampionare e riscrivere dati da file binari
! generati da una simulazione Fortran 77.
! Scritto in Fortran 2008.
!
! -----------------------------------------------------------------------------
! Compilazione (con gfortran):
!   gfortran -o subsampler src/subsample_data.f90 -Wall -Wextra
!
! Esecuzione:
!   ./subsampler <timestep> [<i_stride> <j_stride> <k_stride>]
!
! Esempi:
!   1. Con parametri di default (2, 4, 2):
!      ./subsampler 000100
!      -> Legge 'boum000100.q' e scrive 'sub000100.q'
!
!   2. Con parametri personalizzati (es. 2, 2, 2):
!      ./subsampler 000100 2 2 2
!      -> Legge 'boum000100.q' e scrive 'sub000100.q'
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

  ! Gestione file e argomenti CLI
  character(len=256) :: executable_name
  integer            :: name_length
  character(len=256) :: grid_in_file, data_in_file
  character(len=256) :: grid_out_file, data_out_file
  character(len=10)  :: time_step_str
  character(len=10)  :: stride_str
  integer            :: io_status, dummy_int, arg_count
  integer, parameter :: in_unit = 10, out_unit = 11

  ! Stride per il sottocampionamento
  integer :: i_stride, j_stride, k_stride

  ! ---------------------------------------------------------------------------
  ! Inizio del programma eseguibile
  ! ---------------------------------------------------------------------------

  ! 1. Analisi degli argomenti da linea di comando
  call get_command_argument(0, executable_name, name_length, io_status) ! Ottiene il nome dell'eseguibile
  arg_count = command_argument_count()
  if ((arg_count < 1).or.(arg_count == 2).or.(arg_count == 3).or.(arg_count > 4)) then
    write(error_unit, '("Errore: numero di argomenti non valido.")')
    write(error_unit, '("Uso: ", a, " <timestep_id> [<i_stride> <j_stride> <k_stride>]")') executable_name(:name_length)
    write(error_unit, '("Esempio 1 (default): ", a, " 000100")') executable_name(:name_length)
    write(error_unit, '("Esempio 2 (custom):  ", a, " 000100 2 2 2")') executable_name(:name_length)
    stop 1
  end if

  ! Legge il time step (obbligatorio)
  call get_command_argument(1, time_step_str)

  ! Imposta gli stride: default o da linea di comando
  if (arg_count == 1) then
    i_stride = 2
    j_stride = 4
    k_stride = 2
    write(*, *) "Stride non specificati, si usano i valori di default (2, 4, 2)."
  else ! arg_count == 4
    call get_command_argument(2, stride_str)
    read(stride_str, *) i_stride
    call get_command_argument(3, stride_str)
    read(stride_str, *) j_stride
    call get_command_argument(4, stride_str)
    read(stride_str, *) k_stride
  end if
  write(*, '(a, 3(i0, a))') "Stride utilizzati: (", i_stride, ", ", j_stride, ", ", k_stride, ")"
  write(*, *)

  ! 2. Costruzione dei nomi dei file
  grid_in_file  = 'grid.xyz'
  grid_out_file = 'subgrid.xyz'
  data_in_file  = 'boum' // trim(time_step_str) // '.q'
  data_out_file = 'bous' // trim(time_step_str) // '.q'

  ! ===========================================================================
  ! 3. Processa il file della griglia (grid.xyz)
  ! ===========================================================================
  write(*, '(a, a)') "-> Processando il file della griglia: ", trim(grid_in_file)
  open(unit=in_unit, file=trim(grid_in_file), form='unformatted', status='old', &
       action='read', iostat=io_status)
  if (io_status /= 0) stop "Errore: impossibile aprire il file di input della griglia."
  read(in_unit) dummy_int
  read(in_unit) n1, n2, n3
  write(*, '(a, 3(i0, a))') "   Dimensioni griglia originale: (", n1, ", ", n2, ", ", n3, ")"
  allocate(xxx(n1, n2, n3), yyy(n1, n2, n3), zzz(n1, n2, n3), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per la griglia completa."
  read(in_unit) xxx, yyy, zzz
  close(in_unit)

  ! Esecuzione del sottocampionamento tramite "array slicing"
  n1_new = size(xxx(1::i_stride, 1, 1))
  n2_new = size(xxx(1, 1::j_stride, 1))
  n3_new = size(xxx(1, 1, 1::k_stride))
  write(*, '(a, 3(i0, a))') "   Dimensioni griglia ridotta:   (", n1_new, ", ", n2_new, ", ", n3_new, ")"
  allocate(xxx_new(n1_new, n2_new, n3_new), &
           yyy_new(n1_new, n2_new, n3_new), &
           zzz_new(n1_new, n2_new, n3_new), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per la griglia ridotta."
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
  read(in_unit) dummy_int
  read(in_unit) n1_check, n2_check, n3_check
  ! Verifica che le dimensioni siano consistenti con il file della griglia
  if (n1_check /= n1 .or. n2_check /= n2 .or. n3_check /= n3) stop "Errore: Le dimensioni nel file .q non corrispondono a quelle del file grid.xyz."
  read(in_unit) epsil, epsil, re, sim_time
  allocate(E(n1, n2, n3), vx(n1, n2, n3), vy(n1, n2, n3), &
          vz(n1, n2, n3), pr(n1, n2, n3), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per i dati completi."
  read(in_unit) E, vx, vy, vz, pr
  close(in_unit)

  ! Sottocampionamento dei 5 array
  allocate(E_new(n1_new, n2_new, n3_new),  &
           vx_new(n1_new, n2_new, n3_new), &
           vy_new(n1_new, n2_new, n3_new), &
           vz_new(n1_new, n2_new, n3_new), &
           pr_new(n1_new, n2_new, n3_new), stat=io_status)
  if (io_status /= 0) stop "Errore: impossibile allocare memoria per i dati ridotti."
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
  write(out_unit) epsil, epsil, re, sim_time
  write(out_unit) E_new, vx_new, vy_new, vz_new, pr_new
  close(out_unit)

  deallocate(E_new, vx_new, vy_new, vz_new, pr_new)

  write(*, '(a, a)') "   File dati ridotti creato: ", trim(data_out_file)
  write(*, *)
  write(*, *) "Operazione completata con successo."

end program subsample_data
