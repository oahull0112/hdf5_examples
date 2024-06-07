!
! Number of processes is assumed to be 2
!
PROGRAM hyperslab_by_row_v2

  USE HDF5 ! This module contains all necessary modules 
     
  IMPLICIT NONE

  include 'mpif.h'
  CHARACTER(LEN=10), PARAMETER :: filename = "sds_row.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: dataspace     ! Dataspace identifier in file (likely)
  INTEGER(HID_T) :: plist_id      ! Property list identifier 

  INTEGER :: nbands = 6
  INTEGER ::  ncols = 8
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsf ! Dataset dimensions
                                                    ! in the file.
!     INTEGER, DIMENSION(7) :: dimsfi = (/6,8,0,0,0,0,0/)
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/6,8/)
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsm = (/3,8/) ! Dataset dimensions
                                                    ! in memory.

  INTEGER(HSIZE_T), DIMENSION(2) :: count, count_out
  INTEGER(HSSIZE_T), DIMENSION(2) :: offset, offset_out
  INTEGER(HSIZE_T), DIMENSION(2) :: stride
  INTEGER(HSIZE_T), DIMENSION(2) :: block
  INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
  INTEGER, ALLOCATABLE :: data_out(:,:) ! Data read back in
  INTEGER :: rank = 2 ! Dataset rank 

  INTEGER :: error, error_n  ! Error flags

  ! MPI definitions and calls.
  INTEGER :: mpierror       ! MPI error flag
  INTEGER :: comm, info
  INTEGER :: mpi_size, mpi_rank
  INTEGER :: file_write, ierr, status(MPI_STATUS_SIZE)
  integer, parameter :: send_data_tag = 2001
  integer, parameter :: return_data_tag = 2002

  INTEGER :: i, j

   comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
  CALL MPI_COMM_RANK(comm, mpi_rank, mpierror) 
  ! Quit if mpi_size is not 2

  ! Initialize HDF5 library and Fortran interfaces.
  CALL h5open_f(error) 

  dimsf(1) = nbands
  dimsf(2) = ncols


  call make_file()
  call read_file()


  ! Close FORTRAN interfaces and HDF5 library.
  CALL h5close_f(error)

  CALL MPI_FINALIZE(mpierror)

  contains

    subroutine make_file()

      allocate(data( nbands, ncols) )

      if (mpi_rank .eq. 0) then
        do i = 1, nbands
          do j = 1, ncols
            data(i, j) = i
          end do
        end do

        do i = 1, nbands
            write(*,*) (data(i, j), j = 1, ncols)
          end do

        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
        call h5screate_simple_f(2, dimsf, dataspace, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
          dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsf, error)
        call h5sclose_f(dataspace, error)
        call h5dclose_f(dset_id, error)
        call h5fclose_f(file_id, error)

        file_write = 1

        ! Do an MPI blocking send/receive call to sync up the processes

        write(*,*) "PID: ", mpi_rank, "wrote the file"
        call MPI_send(file_write, 1, MPI_INT, 1, send_data_tag, &
        MPI_COMM_WORLD,ierr)
      
      end if

      if (mpi_rank .eq. 1) then
        call MPI_RECV(file_write, 1, MPI_INT, 0, MPI_ANY_TAG, &
           MPI_COMM_WORLD, status, ierr)
      end if


      
  end subroutine make_file

  subroutine read_file()

    allocate (data_out(dimsm(1), dimsm(2)))
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, mpierror)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, & 
      error, access_prp=plist_id)

    call h5pclose_f(plist_id, error)

    ! bands subroutine stuff starts here
    ! Before adding this subroutine stuff, everything compiled and seemed to be
    ! working, so restart from this point if errors go beyond the point of no
    ! return.


    ! HERE we needd to set the offset, count, etc. Or nothing good will probably
    ! happen.

    count(1) = nbands/mpi_size
    count(2) = 0

    offset(1) = count(1)*mpi_rank ! this may need a -1?
    if (mpi_rank .ne. 0) then
      offset(1) = offset(1) - 1
    end if
   ! offset(1) = 0
    offset(2) = 0

    do j = 1, dimsm(2)
      do i = 1, count(1)
        data_out(i, j ) = mpi_rank
      end do
    end do

    write(*,*) "PID: ", mpi_rank, "has count = ", count(1), &
      "and offset = ", offset(1)

    call h5dopen_f(file_id, dsetname, dset_id, error)

    call h5dget_space_f(dset_id, dataspace, error)
    write(*,*) "dataset id: ", dset_id, "dataspace id: ", dataspace


    call h5screate_simple_f(2, dimsm, memspace, error) ! was dimsm

    count_out(2) = 0
    count_out(1) = count(1)
    offset_out(2) = 0
    offset_out(1) = 0



    ! This two calls were AFTER the hyperslab selections but before the read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error) !
    ! was H5FD_MPIO_INDEPENDENT_F

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, & 
              count, error)

    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_out, &
      count_out, error)

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsm, error, memspace, &
      dataspace, xfer_prp=plist_id) ! this was previously data instead of
    ! data_out

    write(*,*) "PID: ", mpi_rank, " obtained read-in matrix: "

    do i = 1, dimsm(1)
      write(*,*) (data_out(i, j), j=1, dimsm(2))
    end do


    call h5pclose_f(plist_id, error)
    call h5sclose_f(memspace, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dset_id, error)
!
    ! bands subroutine stuff ends here

    call h5fclose_f(file_id, error)

  end subroutine read_file



END PROGRAM hyperslab_by_row_v2
