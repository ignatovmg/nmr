! Copyright (c) 2011-2013 by Evangelos A. Coutsias and Michael J. Wester
! Department of Mathematics and Statistics
! University of New Mexico, Albuquerque, New Mexico, USA
! Written by Michael J. Wester

!==============================================================================
MODULE PDB_coords
!==============================================================================

PUBLIC

integer, PRIVATE, PARAMETER :: DP = kind(1.0d0)

CONTAINS

!==============================================================================

SUBROUTINE get_filenames(filelist, n_files, file)

   implicit none
   character(*), intent(in) :: filelist
   integer, intent(out) :: n_files
   character(*), pointer :: file(:)
   
   character(80) :: line
   integer :: i, k, in, stat;

   in = 1
   open(in, file = filelist, action = 'read', status = 'old', &
            position = 'rewind', iostat = stat)
   if (stat > 0) then
      print '("get_filenames: ",a," cannot be opened.")', trim(filelist)
      stop
   end if

   n_files = 0
   read(in, '(a)', iostat = stat) line
   do while (stat == 0)
      n_files = n_files + 1
      read(in, '(a)', iostat = stat) line
   end do

   close(in)

   allocate(file(n_files))

   open(in, file = filelist, action = 'read', status = 'old', &
            position = 'rewind', iostat = stat)
   do i = 1, n_files
      read(in, '(a)') line

      ! Delete trailing .suffix
      k = index(line, '.', back = .true.) - 1
      file(i) = line(1:k)
   end do
   close(in)

end SUBROUTINE get_filenames

!------------------------------------------------------------------------------

integer FUNCTION n_PDB_coords(PDBfile)

   implicit none
   character(*), intent(in) :: PDBfile

   integer :: in, n, stat
   character(80) :: line

   in = 1
   open(in, file = PDBfile, action = 'read', status = 'old', &
            position = 'rewind', iostat = stat)
   if (stat > 0) then
      print '("n_PDB_coords: ",a," cannot be opened.")', trim(PDBfile)
      stop
   end if

   n = 0
   read(in, '(a)', iostat = stat) line
   do while (stat == 0)
      if (line(1:4) == "ATOM") n = n + 1
      read(in, '(a)', iostat = stat) line
   end do

   close(in)

   n_PDB_coords = n

end FUNCTION n_PDB_coords

!------------------------------------------------------------------------------

SUBROUTINE get_PDB_coords(PDBfile, n_atoms, XYZs)

   implicit none
   character(*), intent(in) :: PDBfile
   integer, intent(in) :: n_atoms
   real (DP), intent(out) :: XYZs(:, :)

   integer :: i, in, stat
   character(80) :: line

   in = 1
   open(in, file = PDBfile, action = 'read', status = 'old', &
            position = 'rewind', iostat = stat)
   if (stat > 0) then
      print '("get_PDB_coords: ",a," cannot be opened.")', trim(PDBfile)
      stop
   end if

   i = 0
   read(in, '(a)', iostat = stat) line
   do while (stat == 0)
      if (line(1:4) == "ATOM") then
         i = i + 1
         if (i <= n_atoms) then
            read(line(31:54), '(3(f8.0))') XYZs(:, i)
         else
            print '("get_PDB_coords: > ",i0," atoms in ",a)', n_atoms, PDBfile
            stop
         end if
      end if
      read(in, '(a)', iostat = stat) line
   end do

   close(in)

end SUBROUTINE get_PDB_coords

! -----------------------------------------------------------------------------

! convert an integer in range 1--9999 into a character of length 4

  character(3) FUNCTION dec3char(ks)
implicit none
  integer ks, k1, k2, k3, k4, kk
      kk = ks
      if (0 .lt. ks .and. 1000 .gt. ks) then
         k1 = (kk - mod(kk, 100))/100
         kk = kk - k1*100
         k2 = (kk - mod(kk, 10))/10
         k3 = kk - k2*10
         k3 = k3+1; k2 = k2+1; k1 = k1+1
         dec3char(1:1) = id(k1); dec3char(2:2) = id(k2)
         dec3char(3:3) = id(k3)
      else
         dec3char(1:1) = 'f'; dec3char(2:2) = 'u';
         dec3char(3:3) = 'z'
      endif
      print*, ks, dec3char(1:3)
  end FUNCTION dec3char

! -----------------------------------------------------------------------------

! Select a one-character identifier to represent the integer k.
character(1) FUNCTION id(k)

   implicit none
   integer, intent(in) :: k

   character(70) :: ids = &
      '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxtz'

   id = ids(k:k)

end FUNCTION id


!==============================================================================
end MODULE PDB_coords
!==============================================================================
