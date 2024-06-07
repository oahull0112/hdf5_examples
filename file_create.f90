PROGRAM FILE_CREATE

  use HDF5
  use mpi

  implicit none

!  include 'mpif.h'
  character(len=10), parameter :: filename = "sds.h5"

  integer(HID_T) :: file_id    ! File identifier
  integer(HID_T) :: plist_id   ! Property list identifier
  integer :: error             ! error flag

  ! MPI definitions and calls:
  integer :: mpierror
  integer :: comm, info
  integer :: mpi_size, mpi_rank
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
  CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)

  ! Initialize hdf5
  CALL h5open_f(error)

  ! Set up file access property list (FAPL) with parallel IO access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)    ! creates FAPL id (plist_id)
  CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)    ! set the comm and info objects to use to open the file later

  ! Create the file collectively
  CALL h5fcreate_f(H5P_FILE_ACCESS_F, plist_id, error)    ! pass in the fapl list to the creation call

  ! Close the property list (FAPL) and the file
  CALL h5pclose_f(plist_id, error)
  CALL h5fclose_f(file_id, error)

  ! Close HDF5 fortran interface
  CALL h5close_f(error)

  ! End MPI
  CALL MPI_FINALIZE(mpierror)

END PROGRAM FILE_CREATE
