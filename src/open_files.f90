

subroutine open_files

  ! takes input file name defined on the command line input, e.g.
  ! %abc.exe run_1.inp grid_file.aba
  ! opens run_1.inp and grid_file.aba plux some additional files
  ! executable name  = command(0)
  ! input file       = command(1)
  ! abaqus grid file = command(2)

  use mod_files
  implicit none
  save

! local variables

  character(len=32), dimension(4) :: command ! command line
  integer                     :: k, km1
  character(len=5), parameter :: out = '.out', aba = '.aba', ver = '.ver'
  logical :: there
  no_arg = command_argument_count()
  write(6,"('no command arguments  = ',i3)") no_arg
  if(no_arg >= 1) then
     call get_command_argument(1,command(1))  ! zeroth argument is program name
     read(command(1),*) inp_file
  end if
  if (no_arg >= 2) then
     call get_command_argument(2,command(2))
     read(command(2),*) aba_file
  end if
  k = index(inp_file,'.')

  !     if file has no extension, then take care of this case

  if (k == 0) k = index(inp_file,' ')
  km1 = k - 1

  !     concatonate non-blank portion of file to create output file

  out_file = inp_file(1:km1) // out          ! output file
  ver_file = inp_file(1:km1) // ver          ! verification file
  if(no_arg >= 1) then
     write(6,"('input file            = ',a32)") inp_file
     open(unit=kinp, file=inp_file, status="old", form='formatted', &
          iostat=ierror) ! input file
     write(6,"('output file           = ',a32)") out_file
     open(unit=kout, file=out_file, status="replace", form='formatted', &
          iostat=ierror) ! output file
  end if

  if (no_arg >= 2) then
     write(6,"('abaqus grid file      = ',a32)") aba_file
     open(unit=kaba, file=aba_file, status="old", form='formatted', &
          iostat=ierror) ! abaqus grid file
  end if

end subroutine open_files
