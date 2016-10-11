!Copyright (c) 2013-2015 by Evangelos A. Coutsias and Michael J. Wester
! Written by EAC 20-8-16

PROGRAM NMR

   USE PDB_coords

   implicit none

   integer, PARAMETER :: DP = kind(1.0d0)

   character(80) :: prefix, filename, filename1, filename2, suffix !,options
   character(80), pointer :: xyz_file(:)
   real (DP), pointer :: XYZs(:, :)
   character(1), pointer :: AtomType(:)
   character(2), pointer :: BondType(:)
   character(3) :: dec3car
   integer done
   integer, pointer :: Bonds(:, :)
   integer, pointer :: ind1(:), ind2(:)
   integer :: i, j, n_atoms, n_bonds, n_files
   real  :: dist, d_min, d_max
   real  :: x1, y1, z1, x2, y2, z2
   integer :: ii, i1, i2, stat, n1, n2

   print*,'Enter Path'
   read '(a)', prefix
!  read '(a)', options

   done = 0
!  
   ii = 0
   write(filename1, '(a,"/files")') trim(prefix)
do while (done .ne. 1)
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

   print*,' enter size of first group of atom indices'
      read*, n1
      allocate(ind1(n1))
   print*,' enter first group of atom indices'
      do i = 1,n1
         print*, 'enter index of 1s grp atom #',i,' of ',n1
         read*, ind1(i)
      end do
   print*,' enter size of second group of atom indices'
      read*, n2
      allocate(ind2(n2))
   print*,' enter first group of atom indices'
      do i = 1,n2
         print*, 'enter index of 2nd grp atom #',i,' of ',n2
         read*, ind2(i)
      end do
   print*,'enter minimum & maximum distance allowed'
   read*, d_min, d_max

   allocate(XYZs(3, n_atoms))

   do j = 1, n_files
      write(filename, '(a,".pdb")'), trim(xyz_file(j))
    ! write(filename, '(a,"/",a,".pdb")') trim(prefix), trim(file(j))
      print*, filename
      call get_PDB_coords(filename, n_atoms, XYZs(:, :))
      x1 = 0.0;y1 = 0.0;z1 = 0.0
      do i = 1,n1
         x1 = x1 + XYZs(1,ind1(i))
         y1 = y1 + XYZs(2,ind1(i))
         z1 = z1 + XYZs(3,ind1(i))
      enddo
      print*,'  x1  ', x1
         x1 = x1/n1;y1 = y1/n1;z1 = z1/n1
      print*,'  x1  ', x1
      x2 = 0.0;y2 = 0.0;z2 = 0.0
      do i = 1,n2
         x2 = x2 + XYZs(1,ind2(i))
         y2 = y2 + XYZs(2,ind2(i))
         z2 = z2 + XYZs(3,ind2(i))
      end do
      print*,'  x2  ', x2
         x2 = x2/n2;y2 = y2/n2;z2 = z2/n2
      print*,'  x2  ', x2
      dist = sqrt((x1 - x2)**2  +  (y1 - y2)**2  +  (z1 - z2)**2)
!     print*, XYZs(1,i1),XYZs(1,i2)
!     print*, XYZs(2,i1),XYZs(2,i2)
!     print*, XYZs(3,i1),XYZs(3,i2)
!     dist = sqrt((XYZs(1,i1) - XYZs(1,i2))**2  +    &
!                 (XYZs(2,i1) - XYZs(2,i2))**2  +    &
!                 (XYZs(3,i1) - XYZs(3,i2))**2)
      print*, dist
      if ((dist .le. d_max) .and. (dist .ge. d_min))  then
         print*, j
         write(9,'(a)',iostat = stat) filename
      endif
   end do
  close(9)
  deallocate(ind1,ind2)
  print*, 'filename1', filename1
  print*, 'filename2', filename2
  print*, 'Done? (1 = yes; 0 = no; -1 repeat)'
  read*, done
! write(filename1, '(a,"/",a")') trim(prefix),filename2
  if (done == 0) then
     filename1 = filename2
  endif
  print*, 'filename1', filename1

end do

end PROGRAM NMR
