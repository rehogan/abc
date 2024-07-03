subroutine read_encl_rad_bc

! read enclosure radiation information

  use mod_constants
  use mod_files
  use mod_global_parameters
  use mod_parsing
  use mod_to_lower
  use mod_rad_encl
  implicit none
  interface
     function compare_string(string,substring)
       implicit none
       logical :: compare_string
       character (len=*), intent (in) :: string, substring
     end function compare_string
  end interface

  character(len=12) :: char1, char2, char3, char4   ! dummy character variables
  character(len=2)  :: delim        ! delimiter
  integer           :: i, j

  rad_encl = 'yes'
  do     ! till you find the record 'end encl rad bc block'
     read(kinp,"(a)") record ! read input record
     write(kout,*) record
     record = to_lower(record)         ! convert to lower case
     if(compare_string(record,'end encl rad')) then
        exit
     else if(compare_string(record,'room width')) then    ! read room width
        read(record,*) char1,char2,delim,rm_w                     
     else if(compare_string(record,'room height')) then   ! read room height
        read(record,*) char1,char2,delim,rm_h                     
     else if(compare_string(record,'nn')) then         ! read encl nodes
        read(record,*) char1,delim, (nn_1(i),i=1,4)
     else if(compare_string(record,'floor')) then      ! read floor
        read(record,*) char1,eps(1)
     else if(compare_string(record,'right wall')) then ! read rt wall
        read(record,*) char1,char2,eps(2)
     else if(compare_string(record,'ceiling')) then    ! read ceiling
        read(record,*) char1,eps(3)
     else if(compare_string(record,'left wall')) then  ! read left wall
        read(record,*) char1,char2,eps(4)
     end if
  end do

  nn_1(5:8) = nn_1(1:4)  ! extend enclosure nodal array
  ! generate the four node sequence for each side
  ! side 1: i  j  k  l
  ! side 2: j  k  l  i
  ! side 3: k  l  i  j
  ! side 4: l  i  j  k
  ! nn_encl_rad(m,n), m is side index
  nn_encl_rad(1,1:4) = nn_1(1:4)
  nn_encl_rad(2,1:4) = nn_1(2:5)
  nn_encl_rad(3,1:4) = nn_1(3:6)
  nn_encl_rad(4,1:4) = nn_1(4:7)

  write(6,"(2es13.5)") rm_w, rm_h
  write(6,"('  eps  nn_i  nn_j  nn_k  nn_l')")
  do i=1,4
     write(6,"(f5.2,4i5)") eps(i),(nn_encl_rad(i,j),j=1,4)
  end do

end subroutine read_encl_rad_bc
