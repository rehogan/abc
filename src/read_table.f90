
subroutine read_table(kfile,max_no_ent,no_ent,tx,ty1,ty2,ty3,ty4,ty5)

! read a table of x-y pairs until "end ..." is encountered
! the independent variable is x
! there can be up to 4 numerical fields plus a charater field
! used for reading temperature dependent thermal properties or
! time dependent boundary conditions.
! the 1st column will be temperature, time, pressure, ...
! while 2nd, 3rd, ... columns are property values or boundary condition values

!    tx(1)          ty1(1)        ty2(1)        ...
!    tx(2)          ty1(2)        ty2(2)        ...
!    ...            ...           ...           ...
!    tx(no_ent)    ty1(no_ent)  ty2(no_ent)  ...

! input:
! kfile      = input file
! max_no_ent = maximum number of table entries

! output:
! no_ent     = number of table entries for this table
! tx         = tabular x (temperature or time) array for this table
! ty1        = 1st property value or boundary condition
! ty2        = 2nd property value or boundary condition
! ty3        = 3rd property value or boundary condition
! ty4        = 4th property value or boundary condition
! ty5        = 5th table entry, this one is a character vector

  use mod_files
  use mod_precision
  use mod_parsing
  implicit none

  interface
     subroutine read_record(kifile,kofile,ncol,ncol_max,record,rec_typ)
       use mod_interface_definitions
       use mod_to_lower
       implicit none
       integer,                 intent(in)  :: kifile, kofile, ncol, ncol_max
       character(len=1),        intent(out) :: rec_typ
       character(len=ncol_max), intent(out) :: record
     end subroutine read_record
  end interface

  real(dp), intent(out), dimension(*) :: tx, ty1
  real(dp), intent(out), dimension(*), optional :: ty2, ty3, ty4
  integer, intent(in)  :: kfile, max_no_ent
  integer, intent(out) :: no_ent
  character(len=1), intent(out), dimension(*), optional :: ty5
  ! local
  integer :: i

  no_ent = 0                       ! initialize data entry counter
  do i=1,max_no_ent                ! loop on number of table entries
     call read_record(kfile,kout,ncol,ncol_max,record,rec_typ)
     if(rec_typ == 'd') then        ! we have valid data record
        if(present(ty5)) then
           read(kfile,*) tx(i),ty1(i),ty2(i),ty3(i),ty4(i),ty5(i)
        else if(present(ty4)) then
           read(kfile,*) tx(i),ty1(i),ty2(i),ty3(i),ty4(i)
        else if(present(ty3)) then
           read(kfile,*) tx(i),ty1(i),ty2(i),ty3(i)
        else if(present(ty2)) then
           read(kfile,*) tx(i),ty1(i),ty2(i)
        else
           read(kfile,*) tx(i),ty1(i)
        end if
        no_ent = no_ent + 1       ! backspace in read_record
     else if(rec_typ == 'e') then   ! end record so exit
        exit
     end if
  end do

  return

end subroutine read_table
