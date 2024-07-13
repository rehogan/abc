
subroutine read_therm_prop

  use mod_constants
  use mod_nodes_elements
  use mod_files
  use mod_mat_prop
  use mod_parsing
!  use mod_time_param
  implicit none
  integer :: i, j
  character(len=12) :: char1, char2, char3   ! dummy character variables
  character(len=2) :: delim        ! delimiter
!  logical :: compare_string
!  external compare_string

  interface
     function compare_string(string,substring)
       implicit none
       logical :: compare_string
       character (len=*), intent (in) :: string, substring
     end function compare_string

     subroutine read_record(kifile,kofile,ncol,ncol_max,record,rec_typ)
       use mod_files
       implicit none
       integer,                 intent(in)  :: kifile, kofile, ncol, ncol_max
       character(len=1),        intent(out) :: rec_typ
       character(len=ncol_max), intent(out) :: record
     end subroutine read_record

     subroutine read_table(kfile,max_no_ent,no_ent,tx,ty1,ty2,ty3,ty4,ty5)
       use mod_files
       use mod_precision
       use mod_parsing
       implicit none
       real(dp), intent(out), dimension(*) :: tx, ty1
       real(dp), intent(out), dimension(*), optional :: ty2, ty3, ty4
       integer, intent(in) :: kfile, max_no_ent
       integer, intent(out) :: no_ent
       character(len=1), intent(out), dimension(*), optional :: ty5
     end subroutine read_table

     subroutine find_string(record,substring)
       use mod_files
       implicit none
       character (len=*), intent (in) :: substring
       character (len=*), intent (out) :: record
       logical :: compare_string
       external compare_string
     end subroutine find_string

  end interface
  ! allocate arrays related to thermal properties
  allocate(no_cp_ent(max_no_mat),   no_cond_ent(max_no_mat), &
           no_emit_ent(max_no_mat), no_abs_ent(max_no_mat))
  allocate(rho(max_no_mat), mat_name(max_no_mat))
  allocate(T_cp(max_no_ent,max_no_mat),     TT_cp(max_no_ent,max_no_mat))
  allocate(T_cond(max_no_ent,max_no_mat),   TT_cond(max_no_ent,max_no_mat))
  allocate(T_emit(max_no_ent,max_no_mat),   TT_emit(max_no_ent,max_no_mat))
  allocate(T_abs(max_no_ent,max_no_mat),    TT_abs(max_no_ent,max_no_mat))
  allocate(T_int_en(max_no_ent,max_no_mat))
  allocate(mat_num(max_no_mat))
! we have already found the 'begin thermal property' record before entering

  do j=1,max_no_mat               ! loop on (unknown) number of materials

     mat_num(j) = j                ! material table numbered in order read in

! look for beginning of mat prop data; be careful to exit at bottom of loop

     call find_string(record,'begin material')
     read(record,*)char1,char2,delim,mat_name(j) ! read material name

     do   ! loop on # material properties, e.g. density, reference enthalpy
        ! conductivity, heat capacity, emittance, absorptance, in any order
        call read_record(kinp,kout,ncol,ncol_max,record,rec_typ)
        if(rec_typ == 'e') exit         ! end of this material block
        if(compare_string(record,'begin density')) then      ! density
           call read_record(kinp,kout,ncol,ncol_max,record,rec_typ) ! find data record
           read(kinp,*)char1,delim,rho(j)       ! read density value
           call find_string(record,'end density') ! locate end density record
        else if(compare_string(record,'begin heat')) then      ! heat capacity
           call read_table(kinp,max_no_ent,no_cp_ent(j),tt_cp(1,j),t_cp(1,j))
        else if(compare_string(record,'begin cond')) then ! conductivity
           call read_table(kinp,max_no_ent,no_cond_ent(j),tt_cond(1,j),t_cond(1,j))
        else if(compare_string(record,'begin emit')) then ! emittance
           call read_table(kinp,max_no_ent,no_emit_ent(j),tt_emit(1,j),t_emit(1,j))
        else if(compare_string(record,'begin abs')) then  ! absorptance
           call read_table(kinp,max_no_ent,no_abs_ent(j),tt_abs(1,j),t_abs(1,j))
        else
           write(6,"('****invalid material property type')")
           stop
        end if

     end do   ! end of loop on # material properties for material j

! test to see if we are "end material" or "end thermal property block"
! if not, then we will read another material

     call read_record(kinp,kout,ncol,ncol_max,record,rec_typ)

     if(compare_string(record,'thermal property')) then
        exit
     else if(compare_string(record,'begin material')) then
        backspace(kinp)
     end if

! decide if you want to force alpha = eps in radiative property tables
! based on first temperature entry in table for material j
! t_abs and tt_abs were initialized to zero in subroutine defaults

     if(abs(tt_abs(1,j)) < v_small) then
        do i=1,no_emit_ent(j)
           t_abs(i,j) = t_emit(i,j)    ! abs = emit
           tt_abs(i,j) = tt_emit(i,j)  ! absorbtance t = emittance t
        end do
     no_abs_ent(j) = no_emit_ent(j)
     end if

  end do    ! end of looop on # materials

  num_mat = j               ! set number of materials

  return

end subroutine read_therm_prop
