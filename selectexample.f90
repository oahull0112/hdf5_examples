PROGRAM SELECTEXAMPLE

  USE HDF5 ! This module contains all necessary modules 
         
  IMPLICIT NONE
  ! RE-WORK SO IT RUNS IN MPI PARALLEL! (major)
  ! IF want to eventually make this a module that can do the timings, just make this a subroutine that takes nbands, ncols, and
  ! excl_bands as input.
 
  CHARACTER(LEN=7), PARAMETER :: filename = "sdsf.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name
  CHARACTER(LEN=7), PARAMETER :: file_out = "fout.h5"  ! File name of output after selections
  CHARACTER(LEN=6), PARAMETER :: dsetname_out = "IntOut"
 
  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: fout_id       ! Output file id
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
  INTEGER(HID_T) :: memspace      ! memspace identifier 
 
  INTEGER :: nbands = 10
  INTEGER :: ncols  = 7 
 
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsf ! Dataset dimensions. ** This will be defined in "main"
  INTEGER(HSIZE_T), DIMENSION(2) :: count ! Size of the hyperslab in the file
  INTEGER(HSIZE_T), DIMENSION(2) :: offset ! Hyperslab offset in the file
  INTEGER(HSIZE_T), DIMENSION(2) :: count_out ! Size of the hyperslab in memory
  INTEGER(HSIZE_T), DIMENSION(2) :: offset_out = (/0,0/) ! hyperslab offset in memory (first row you want to keep has no offset)
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsm ! = (/6, 6/) ! (/10,6/) ! Dataset dimensions in memory ** MUST CHANGE
  INTEGER(HSIZE_T), DIMENSION(2) :: dims_out ! Buffer to read in dataset dimensions
  INTEGER, DIMENSION(1) :: nbands_excl_temp
  INTEGER :: nbands_excl, ncore_excl
 
  INTEGER, ALLOCATABLE :: data(:,:)
  INTEGER, ALLOCATABLE :: data_out(:,:)
  INTEGER :: dsetrank = 2 ! Dataset rank ( in file )
  INTEGER :: memrank = 2  ! Dataset rank ( in memory )
  INTEGER :: i, j, k, q
 
  INTEGER :: error  ! Error flag
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
 
  REAL :: start, finish, diff, innerstart, innerfinish, innerdiff ! timings
 
  INTEGER, DIMENSION(4) :: excl_bands ! ** Want to take this as input eventually
  call cpu_time(start)
  excl_bands(1) = 1
 !     do q = 2, 500
 !       excl_bands(q) = q + 1
 !     enddo
 ! uncomment to generate bands automatically, or define manually as in next
 ! line: 
  excl_bands = (/2,3, 7,11/)
 
  allocate( data(nbands, ncols) )
 
  dimsf(1) = nbands
  dimsf(2) = ncols
 
  ! These won't change. Just setting the dimension we don't care about
  count(2) = ncols
  count_out(2) = ncols
  offset(2) = 0
 
  nbands_excl_temp = SHAPE(excl_bands)
  nbands_excl = nbands_excl_temp(1)
  ncore_excl = 0
  
  ! the dimensions of what you want to read to memory -- based on dimsf and excluded bands here
  dimsm(2) = ncols
  dimsm(1) = nbands - nbands_excl + 1 ! + 1 because excl_bands has the last term as nbands, which isn't actually excluded
  allocate( data_out(dimsm(1), dimsm(2)) )
 
  ! dimensions of data you're writing to file (matches with build_data subroutine if you want the whole matrix written to file)
  data_dims(1) = nbands
  data_dims(2) = ncols
 
  ! Determine number of excluded core states:
  call determine_ncore_excl()
 
  print *, ncore_excl
 
 ! Write data to the HDF5 file.  
 
  ! Make some matrix, then write it to an HDF5 file
  ! Modify the build_data subroutine if you want the matrix elements to look different.
  call build_data()
 
  write(*,*) " "
  ! Initialize FORTRAN interface. 
  CALL h5open_f(error) 
  ! Create a new file using default properties.
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
   ! created, into a 2-dimensional array.
 
     ! Open the file
  CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
      ! Open the  dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, error)
  ! Get dataset's dataspace identifier.
  CALL h5dget_space_f(dset_id, dataspace, error)
  call cpu_time(innerstart)
  ! Select hyperslab in the dataset.

  ! This is the "main" subroutine. Generates the matrix that includes only the
  ! bands you wanted to include, and reads it in to data_out
  call exclude_bands()
 
  ! Display the new matrix
  do i = 1, dimsm(1)
      print *, (data_out(i,j), j = 1, dimsm(2))
  end do
 
  ! Close the dataspace for the dataset.
  CALL h5sclose_f(dataspace, error)
  ! Close the memoryspace.
  CALL h5sclose_f(memspace, error)
  ! Close the dataset.
  CALL h5dclose_f(dset_id, error)
  ! Close the file.
  CALL h5fclose_f(file_id, error)
 
  ! Write the new matrix (data_out) into a new file.
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

      if (excl_bands(1) .eq. 1) then ! if first excluded band is the first overall band, then have at least 1 excluded core state
        ncore_excl = ncore_excl + 1 
        do i = 1, nbands_excl-1
          if (excl_bands(i + 1) - excl_bands(i) .ne. 1) exit
          ncore_excl = ncore_excl + 1
        enddo
      endif

    end subroutine determine_ncore_excl

    subroutine build_data()

      do i = 1, nbands
           do j = 1, ncols
               data(i,j) = i
                !data(i,j) = 1 + (i-1) + (j-1);
           end do
      end do

      do i = 1, nbands
          print *, (data(i,j), j = 1, ncols)
      end do

    end subroutine build_data

    subroutine initialize_data_out()

     ! Initialize data_out array.
      do j = 1, dimsm(2) ! dimensions ought to match the memory space dimensions
        do i = 1, dimsm(1) 
              data_out(i,j) = 0;
          end do
      end do

    end subroutine initialize_data_out

    subroutine exclude_bands()

      offset(1) = ncore_excl
      if (ncore_excl .eq. 0) then
        count(1) = excl_bands(1) - 1
      else
        count(1) = excl_bands(ncore_excl + 1) - excl_bands(ncore_excl) - 1
      end if
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error) 
      
      do i = 1, ( nbands_excl - ncore_excl) - 1
        offset(1) = excl_bands(ncore_excl + i) !+ 1
        count(1) = excl_bands(ncore_excl + i + 1) - excl_bands(ncore_excl + i) - 1
        CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_OR_F, offset, count, error)
      enddo  
 
      ! Create memory dataspace.
      CALL h5screate_simple_f(memrank, dimsm, memspace, error)
      ! Select hyperslab in memory.
      
      if (ncore_excl .eq. 0) then
        count_out(1) = excl_bands(1) - 1
      else
        count_out(1) = excl_bands(ncore_excl + 1) - excl_bands(ncore_excl) - 1
      end if
      ! offset_out does not change for the first call!
      CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_out, count_out, error) 
      
      offset_out(1) = count_out(1)
      do i = 1, (nbands_excl - ncore_excl) - 1
        count_out(1) = excl_bands(ncore_excl + i + 1) - excl_bands(ncore_excl + i) - 1
        CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_OR_F, offset_out, count_out, error)
        offset_out(1) = offset_out(1) + count_out(1)
      enddo
 
      call cpu_time(innerfinish)
 
      ! Read data from hyperslab in the file into the hyperslab in memory
      ! Note it only reads in the stuff you wanted, as specified in the
      ! hyperslab union
      data_dims(1) = dimsm(1) 
      data_dims(2) = dimsm(2) 
      CALL H5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, data_dims, error, memspace, dataspace)
      ! Display data_out array
 
    end subroutine exclude_bands

     END PROGRAM SELECTEXAMPLE 
