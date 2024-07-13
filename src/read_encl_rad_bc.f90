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
     else if(compare_string(record,'nn')) then            ! read encl nodes
        read(record,*) char1,delim, (nn_encl_rad(i),i=1,4)
     else if(compare_string(record,'floor')) then         ! read floor
        read(record,*) char1,eps(1)
     else if(compare_string(record,'right wall')) then    ! read rt wall
        read(record,*) char1,char2,eps(2)
     else if(compare_string(record,'ceiling')) then       ! read ceiling
        read(record,*) char1,eps(3)
     else if(compare_string(record,'left wall')) then     ! read left wall
        read(record,*) char1,char2,eps(4)
     end if
  end do

  ! generate the two node sequence for each side
  ! side 1: i  j  
  ! side 2: j  k  
  ! side 3: k  l  
  ! side 4: l  i  

  surf_ij(1,1:2) = nn_encl_rad(1:2) ! surf 1
  surf_ij(2,1:2) = nn_encl_rad(2:3) ! surf 2
  surf_ij(3,1:2) = nn_encl_rad(3:4) ! surf 3
  surf_ij(4,1)   = nn_encl_rad(4)   ! surf 4, node 1
  surf_ij(4,2)   = nn_encl_rad(1)   ! surf 4, node 2

  write(6,"(2es13.5)") rm_w, rm_h
  write(6,"('  eps   side   nn_i  nn_j ')")
  do i=1,4
     write(6,"(f5.2,3i5)") eps(i),i,(surf_ij(i,j),j=1,2)
  end do

end subroutine read_encl_rad_bc
