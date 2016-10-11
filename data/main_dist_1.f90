!Copyright (c) 2013-2015 by Evangelos A. Coutsias and Michael J. Wester
! Written by EAC 20-8-16

PROGRAM NMR

   USE PDB_coords

   implicit none

   integer, PARAMETER :: DP = kind(1.0d0)

   character(80) :: prefix, filename, filename1, filename2, suffix !,options
   character(80), pointer :: xyz_file(:)
   real (DP), pointer :: XYZs(:, :, :)
   character(1), pointer :: AtomType(:)
   character(2), pointer :: BondType(:)
   character(3) :: dec3car
   integer done
   integer, pointer :: Bonds(:, :)
   integer :: i, n_atoms, n_bonds, n_files
   real  :: dist, d_min, d_max
   integer :: ii, i1, i2, stat

   print*,'Enter Path'
   read '(a)', prefix
!  read '(a)', options

   done = 0
!  
   ii = 0
   write(filename1, '(a,"/files")') trim(prefix)
do while (done == 0)
   ii = ii + 1
   suffix = dec3char(ii)
   print*, suffix
   write(filename2, '("files",a)') trim(suffix)
   print*, filename2
   open(9, file = filename2, action = 'write', iostat = stat)
!  write(filename, '(a,"/files")') trim(prefix)
   call get_filenames(filename1, n_files, xyz_file)
   print '("n_files = ",i0)', n_files
   print*, xyz_file

  !write(filename, '(a,"/",a)') trim(prefix), xyz_file(1)
  !print*, xyz_file(1)

   write(filename, '(a,"/",a,".pdb")') trim(prefix), trim(xyz_file(1))
   print*, filename
   n_atoms = n_PDB_coords(filename)
   print '("n_atoms = ",i0)', n_atoms

   print*,' enter pair of atom indices'
   read*, i1, i2
   print*,'enter minimum & maximum distance allowed'
   read*, d_min, d_max

   allocate(XYZs(3, n_atoms, n_files))

   do i = 1, n_files
      write(filename, '(a,".pdb")'), trim(xyz_file(i))
    ! write(filename, '(a,"/",a,".pdb")') trim(prefix), trim(file(i))
      print*, filename
      call get_PDB_coords(filename, n_atoms, XYZs(:, :, i))
      print*, XYZs(1,i1,i),XYZs(1,i2,i)
      print*, XYZs(2,i1,i),XYZs(2,i2,i)
      print*, XYZs(3,i1,i),XYZs(3,i2,i)
      dist = sqrt((XYZs(1,i1,i) - XYZs(1,i2,i))**2  +    &
                  (XYZs(2,i1,i) - XYZs(2,i2,i))**2  +    &
                  (XYZs(3,i1,i) - XYZs(3,i2,i))**2)
      print*, dist
      if ((dist .le. d_max) .and. (dist .ge. d_min))  then
         print*, i
         write(9,'(a)',iostat = stat) filename
      endif
   end do
  close(9)
  print*, 'filename1', filename1
  print*, 'filename2', filename2
! write(filename1, '(a,"/",a")') trim(prefix),filename2
  filename1 = filename2
  print*, 'filename1', filename1
  print*, 'Done? (1 = yes; 0 = no)'
  read*, done

end do

end PROGRAM NMR
