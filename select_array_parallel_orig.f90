PROGRAM select_array_parallel

  USE HDF5 ! This module contains all necessary modules 
         
  IMPLICIT NONE

  include 'mpif.h'
  ! RE-WORK SO IT RUNS IN MPI PARALLEL! (major)
  ! IF want to eventually make this a module that can do the timings, just make this a subroutine that takes nbands, ncols, and
  ! excl_bands as input.
  ! Specify band ranges to include in incl_bands

  CHARACTER(LEN=7), PARAMETER :: filename = "sdsf.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name
  CHARACTER(LEN=7), PARAMETER :: file_out = "fout.h5"  ! File name of output after selections
  CHARACTER(LEN=6), PARAMETER :: dsetname_out = "IntOut"
 
  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: fout_id       ! Output file id
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
  INTEGER(HID_T) :: memspace      ! memspace identifier 
  INTEGER(HID_T) :: plist_id      ! Property list identifier
 
  INTEGER :: nbands = 10000
  INTEGER :: ncols  = 10000 ! Can change this later for the testing (make this dimension very large, measure timings)
 
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsf ! Dataset dimensions. ** This will be defined in "main"
  INTEGER(HSIZE_T), DIMENSION(2) :: count ! Size of the hyperslab in the file
  INTEGER(HSIZE_T), DIMENSION(2) :: offset ! Hyperslab offset in the file
  INTEGER(HSIZE_T), DIMENSION(2) :: count_out ! Size of the hyperslab in memory
  INTEGER(HSIZE_T), DIMENSION(2) :: offset_out = (/0,0/) ! hyperslab offset in memory (first row you want to keep has no offset)
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsm ! = (/6, 6/) ! (/10,6/) ! Dataset dimensions in memory ** MUST CHANGE
  INTEGER(HSIZE_T), DIMENSION(2) :: dims_out ! Buffer to read in dataset dimensions
!  INTEGER, DIMENSION(1) :: nbands_excl_temp
 
  INTEGER, ALLOCATABLE :: data(:,:)
  INTEGER, ALLOCATABLE :: data_out(:,:)
  INTEGER :: dsetrank = 2 ! Dataset rank ( in file )
  INTEGER :: memrank = 2  ! Dataset rank ( in memory )
  INTEGER :: i, j, k, q
 
  INTEGER :: error  ! Error flag
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  INTEGER, DIMENSION(2) :: incl_size ! nbands_incl(2) = 2 always. Care about nbands_incl(1)
!  INTEGER, DIMENSION(4, 2) :: incl_bands
  INTEGER, ALLOCATABLE :: incl_bands(:, :)

  INTEGER :: nbands_excl, ncore_excl, nbands_incl

  REAL :: start, finish, diff, innerstart, innerfinish, innerdiff ! timings
 
  INTEGER, DIMENSION(4) :: excl_bands ! ** Want to take this as input eventually


  INTEGER :: mpierror, comm, info,mpi_size, mpi_rank

  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  call MPI_INIT(mpierror)
  call MPI_COMM_SIZE(comm, mpi_size, mpierror)
  call MPI_COMM_RANK(comm, mpi_rank, mpierror)

  call cpu_time(start)
 
!  allocate( incl_bands( (nbands/2) - 1, 2) )
  allocate (incl_bands( 1, 2))
  incl_bands(1, 1) = 5001
  incl_bands(1, 2) = 10000


!  incl_bands(1, 1) = 2
!  incl_bands(1, 2) = 2
!  incl_bands(2, 1) = 4
!  incl_bands(2, 2) = 6
!  incl_bands(3, 1) = 8
!  incl_bands(3, 2) = 9
!  incl_bands(4, 1) = 20
!  incl_bands(4, 2) = 25 

!  do i = 1, nbands/2 - 1
!    incl_bands(i, 1) = 2*i
!    incl_bands(i, 2) = 2*i
!  end do


 
  incl_size = shape(incl_bands)
 
  allocate( data(nbands, ncols) )
 
  dimsf(1) = nbands
  dimsf(2) = ncols
 
  ! These won't change. Just setting the dimension we don't care about
  count(2) = ncols
  count_out(2) = ncols
  offset(2) = 0
 
  
  ! the dimensions of what you want to read to memory -- based on dimsf and excluded bands here
 
  call determine_ncore_excl()
  call determine_nbands_excl()
  call determine_nbands_incl()

  write(*,*) ncore_excl, nbands_excl, nbands_incl


   dimsm(2) = ncols
   dimsm(1) = nbands_incl 
   print *, ncore_excl

 
  allocate( data_out(dimsm(1), dimsm(2)) )
 
  ! dimensions of data you're writing to file (matches with build_data subroutine if you want the whole matrix written to file)
!  data_dims(1) = nbands
!  data_dims(2) = ncols
 ! Write data to the HDF5 file.  
 
  ! Make some matrix, then write it to an HDF5 file
  call build_data()
 
  write(*,*) " "
  ! Initialize FORTRAN interface. 
  CALL h5open_f(error) 

  ! Set up file access property list with parallel I/O access
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)


  ! Create a new file using default properties.

  ! Don't actually want to create a parallel file here, just want to read one
  ! in...
  ! MONDAY: follow along between BGW code and the hyperslab_by_row sample code
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
  ! Create the data space for the  dataset. 
  CALL h5screate_simple_f(dsetrank, dimsf, dataspace, error)
  ! Create the dataset with default properties.
  CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, dset_id, error)
  ! Write the dataset.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, data_dims, error)
  ! Close the dataspace for the dataset.
  CALL h5sclose_f(dataspace, error)
  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  ! Close the file.
  CALL h5fclose_f(file_id, error)
 
   ! This  part of the code reads the hyperslab from the sds.h5 file just 
   ! created, into a 2-dimensional plane of the 3-dimensional dataset.
 
     ! Open the file
  CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
      ! Open the  dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, error)
  ! Get dataset's dataspace identifier.
  CALL h5dget_space_f(dset_id, dataspace, error)
  call cpu_time(innerstart)
  ! Select hyperslab in the dataset.

  call exclude_bands()
 
  ! Display the new matrix
 ! do i = 1, dimsm(1) ! ** CHANGE
 !     print *, (data_out(i,j), j = 1, dimsm(2)) ! ** CHANGE
 ! end do
 
  ! Close the dataspace for the dataset.
  CALL h5sclose_f(dataspace, error)
  ! Close the memoryspace.
  CALL h5sclose_f(memspace, error)
  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  ! Close the file.
  CALL h5fclose_f(file_id, error)
 
  ! WRITE THE SELECTED DATA TO A NEW FILE:
  CALL h5fcreate_f(file_out, H5F_ACC_TRUNC_F, file_id, error)
  CALL h5screate_simple_f(dsetrank, dimsm, dataspace, error)
  CALL h5dcreate_f(file_id, dsetname_out, H5T_NATIVE_INTEGER, dataspace, dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_out, data_dims, error)
  call h5sclose_f(dataspace, error)
  call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id, error)
 
  ! Close FORTRAN interface.
  CALL h5close_f(error)
 
  call cpu_time(finish)
  innerdiff = innerfinish - innerstart
  diff = finish - start
 
  print *, "Total time: ", diff
  print *, "Hyperslabbing time: ", innerdiff
 
  deallocate(data_out)
  deallocate(data)
 
  contains

    subroutine determine_ncore_excl()

      ncore_excl = incl_bands(1,1) - 1

    end subroutine determine_ncore_excl

    subroutine determine_nbands_excl()

      ! OR, bands_excl = nbands - nbands_incl ...

      nbands_excl = 0
      do i = 1, incl_size(1) - 1
        nbands_excl = nbands_excl + incl_bands(i+1, 1) - incl_bands(i, 2) - 1
      end do

      nbands_excl = nbands_excl + ncore_excl + (nbands - incl_bands( incl_size(1), incl_size(2) ))

    end subroutine determine_nbands_excl

    subroutine determine_nbands_incl()

      nbands_incl = 0
      do i = 1, incl_size(1)
        nbands_incl = nbands_incl + incl_bands(i, 2) - incl_bands(i, 1) + 1
      end do

      ! OR, nbands_incl = nbands - nbands_excl ...

    end subroutine determine_nbands_incl

    subroutine build_data()

      do i = 1, nbands
           do j = 1, ncols
               data(i,j) = i
                !data(i,j) = 1 + (i-1) + (j-1);
           end do
      end do

!      do i = 1, nbands
!          print *, (data(i,j), j = 1, ncols)
!      end do

    end subroutine build_data

    subroutine initialize_data_out()

     ! Initialize data_out array.
      do j = 1, dimsm(2) ! ** MUST CHANGE
        do i = 1, dimsm(1) ! 10 ** MUST CHANGE if this works then dimsf(1) = 6?
              data_out(i,j) = 0
          end do
      end do

    end subroutine initialize_data_out

    subroutine exclude_bands()

      offset(1) = incl_bands(1,1) - 1! see if this should be ncore_excl or ncore_excl + 1
      count(1) = incl_bands(1,2) - incl_bands(1,1) + 1
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error) 
      
      do i = 2, incl_size(1)
        offset(1) = incl_bands(i, 1) - 1
        count(1) = incl_bands(i, 2) - incl_bands(i, 1) + 1 

        CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_OR_F, offset, count, error)
      enddo  
 
      ! Create memory dataspace.
      CALL h5screate_simple_f(memrank, dimsm, memspace, error)
      ! Select hyperslab in memory.
      
      count_out(1) = incl_bands(1,2) - incl_bands(1,1) + 1

      ! offset_out does not change for the first call!
      CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_out, count_out, error) 
      
      offset_out(1) = count_out(1)
      do i = 2, incl_size(1) 
        count_out(1) = incl_bands(i, 2) - incl_bands(i, 1) + 1
        CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_OR_F, offset_out, count_out, error)
        offset_out(1) = offset_out(1) + count_out(1)
      enddo
 
      call cpu_time(innerfinish)
 
      ! Read data from hyperslab in the file into the hyperslab in memory and display.
      data_dims(1) = dimsm(1) 
      data_dims(2) = dimsm(2) 
      CALL H5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, data_dims, error, memspace, dataspace)
      ! Display data_out array
 
    end subroutine exclude_bands

END PROGRAM select_array_parallel
