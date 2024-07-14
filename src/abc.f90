! this is test
! add 2nd line
program abc
  ! abc = any body's code
  use mod_precision
  use mod_constants
  use mod_bc
  use mod_global_parameters
  use mod_nodes_elements
  use mod_files
  use mod_time_params
  use mod_rad_encl
  use mod_rad_encl_lib
  implicit none
  integer                   :: i, j, iter, b_el_no
  integer, dimension(2)     :: nn_b
  real(dp)                  :: del_z, a, h, dhdt, t_inf, exp
  real(dp)                  :: qdot_pp, dqdot_pp_dT
  real(dp), dimension(4,4)  :: U_mat      ! {q_r} = [U]{e_b}
  character, dimension(15)  :: title
  real(dp),  dimension(200) :: x_num, T_num
  integer                   :: nmod
  character (len=11)        :: geom_typ   ! verification only
  character (len=7)         :: geom

  nmod(i) = 1 + i - (i/3)*3  ! function to cycle time step 
                             ! indices 1,2,3,1,2,3,...
  interface
     subroutine air_tp_props(temp, h, dhdt)
       use mod_precision
       use mod_constants
       real(kind=dp), intent(in)  :: temp
       real(kind=dp), intent(out) :: h, dhdt
     end subroutine air_tp_props

     subroutine bc_del_st(node_no)
       ! node_no = node number of specified temperature node
       ! T_spec  = temperature of specified node
       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       integer,        intent(in) :: node_no
     end subroutine bc_del_st

     subroutine band(n,nfb,c,b,x)
       use mod_precision
       implicit none
       integer,                         intent(in)    :: n, nfb    !! nk, nl
       real(kind=dp),  dimension(:,:),  intent(out)   :: c
       real(kind=dp),  dimension(:),    intent(out)   :: b, x
     end subroutine band


     subroutine initialize
       ! no arguments
     end subroutine initialize

     subroutine read_aba_node_elem_blk
       ! no arguments
     end subroutine read_aba_node_elem_blk

     subroutine cond_cap
       ! steady state conduction w/o any sources but with temp dependent cond
       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       integer            :: j, jj, k, kk, mat_no_eb
       integer, parameter :: nnpe = 4
       real(kind=dp),  dimension(4,4)              :: el_cond  ! el cond matrix
     end subroutine cond_cap

     subroutine bc_del_ff_rad(b_el_no,nn_b,eps,fij,t_r)
       ! b_el_no = boundary elem no
       ! i,j     = node no that define side length S_(i,j) for this b_el_no
       ! T_r     = temp to which surf radiates to

       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       integer,       intent(in)               :: b_el_no
       integer,       intent(in), dimension(2) :: nn_b
       real(kind=dp), intent(in)               :: eps, fij, t_r
     end subroutine bc_del_ff_rad

     subroutine bc_del_natl_conv(geom,b_el_no,nn_b,t_inf)
       use mod_precision
       use mod_constants
       use mod_bc
       use mod_nodes_elements
       implicit none
       character,     intent(in), dimension(7) :: geom
       integer,       intent(in)               :: b_el_no
       integer,       intent(in), dimension(2) :: nn_b
       real(kind=dp), intent(in)               :: t_inf
     end subroutine bc_del_natl_conv

     subroutine bc_del_nl_flux(b_el_no,nn_b,a,exp,t_inf)
       ! non-linear flux of the form q_dot^" = a(T_inf - T_w)**m
       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       integer,       intent(in)               :: b_el_no
       integer,       intent(in), dimension(2) :: nn_b
       real(kind=dp), intent(in)               :: a, exp, t_inf
     end subroutine bc_del_nl_flux

     subroutine bc_del_spec_flux(b_el_no,nn_b,del_z,qdot_pp,dqdot_pp_dT)
       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       integer,       intent(in)               :: b_el_no
       integer,       intent(in), dimension(2) :: nn_b
       real(kind=dp), intent(in)               :: del_z, qdot_pp, dqdot_pp_dT
     end subroutine bc_del_spec_flux

!!$     subroutine get_encl_u_mat(U_mat)
!!$       ! {q_r} = [U]{e_b}, [U] = U_mat
!!$       use mod_precision
!!$       use mod_constants
!!$       use mod_rad_encl
!!$       implicit none
!!$       real(dp),  intent(out),  dimension(ns,ns) :: U_mat  ! {q_r} = [U]{e_b}
!!$     end subroutine get_encl_u_mat

     subroutine read_mat_prop_data
       use mod_mat_prop
       implicit none
     end subroutine read_mat_prop_data

     subroutine read_input
       use mod_files
       use mod_parsing
       implicit none
     end subroutine read_input

     subroutine verification(num_node,geom_typ,n_time,time,r,x,z,t,sdot)
       use mod_precision
       use mod_constants
       implicit none
       integer,   intent(in)               :: num_node, n_time
       real(dp),  intent(in), optional     :: time, sdot
       real(dp),  intent(in), dimension(*) :: r, t, x, z
       character(len=11), intent(in)     :: geom_typ
     end subroutine verification

     subroutine open_files
       use mod_files
       implicit none
       save
     end subroutine open_files

!!$     subroutine vfac(w,h,f,s)
!!$       use mod_precision
!!$       use mod_constants
!!$       implicit none                                        
!!$       real(dp), intent(in)                  :: w, h
!!$       real(dp), intent(out), dimension(6)   :: s
!!$       real(dp), intent(out), dimension(4,4) :: f
!!$     end subroutine vfac

  end interface
  
  len_conv = one
  call open_files              ! open input, output, abaqua, ... files
  call read_input              ! read input file
  if(rad_encl == 'yes') then   
     call vfac(rm_w,rm_h,f,s)  ! calculate and save view factors for rm encl
     call get_encl_u_mat(U_mat)! {q_r} = [U][e_b}, U_mat = [U]
  end if     
  if(no_arg <= 1) stop         ! allows testing of input reader without grid data
  call read_aba_node_elem_blk  ! read abaqus data file for node/element data
  call initialize              ! initialize some arrays, do geometry
  time_ini = zero              ! initial time
  time_fin = zero              ! final time
  dtime(:) = zero
  if (abs(time_ini - time_fin) < small) then
     ss = 'yes'
     last_index = first_index
  end if

  do n_time = first_index,last_index
     theta = one                  ! 1st time step must have theta = 1
     nm1 = nmod(n_time - 1)       ! time index n-1
     n = nmod(n_time)             ! time index n
     np1 = nmod(n_time + 1)       ! time index n+1
     time = time + dtime(np1)     ! increment problem time 

     do iter=1,10                 ! begin iteration loop
        call cond_cap             ! conduction + capaticance

        ! loop on no convective boundary conditions
        geom = 'verif'
        do j=1,no_conv_bc
           b_el_no = ne_bc_conv(j)  ! global bc elem no
           nn_b(1) = ni_bc_conv(j)  ! global bc node no i
           nn_b(2) = nj_bc_conv(j)  ! global bc node no j
           t_inf = t_inf_cbc(j)
           call bc_del_natl_conv(geom,b_el_no,nn_b,t_inf)
        end do

        if(rad_encl == 'yes') then   
           call bc_del_encl_rad(U_mat)! {q_r} = [U][e_b}, U_mat = [U]
        end if
        

        ! specified temperature bc's
        do j=1,no_spec_t_bc
           call bc_del_st(ni_bc_spec_t(j))
        end do

        ! solve linear system of equations
        call band(no_nodes,nfb,c,b,z) ! solve linear system
        t_li(:) = t_li(:) + z(:)      ! add corr vec to last iter

        ! calculate maximum temperature correction
        max_corr = zero
        do j=1,no_nodes
           max_corr = max(max_corr,abs(z(j)))
        end do

        ! beginning of verification block
        ! warning: this section is customized for the specific verification
        ! problem. present problem assumes the following numbering system 
        ! for the nodes:
        !   1____________3_________________5__
        !  |             |                 |
        !  |             |                 |
        !  |             |                 |
        !  |_____________|_________________|
        !   2            4                 6

        i = 0
        do j=1,no_nodes,2
           i = i + 1
           T_num(i) = T_li(j)
           x_num(i) = x(j)
        end do
        write(kver,"('# iter = ',i3)") iter
        geom_typ = 'PLANAR'
        if(verif == 'yes') then
           call verification(no_nodes/2,geom_typ,n_time,time,&
                x_num,x_num,x_num,T_li)
        end if
        ! write output
             write(6,"(' time = ',es13.5, '   iter = ',i3)") time, iter
             write(6,"('   i          T_li                DT')") 
             i = 0
             do j=1,no_nodes,2
                i = i + 1
                write(6,"(i4,es23.15,es15.7)") j,t_li(j),z(j)
             end do
        ! zero linear eqn matrix/vector for next iteration
        c(:,:) = zero
        b(:) = zero
        z(:) = zero

        ! convergence check for non-linear iteration
        if(max_corr < v_small) exit
     end do    ! end of non-linear iteration loop
     ! fill temperature arrays in prep for next time step
     T(:,np1) = T_li(:)
     Tdot(:,np1) = (T(:,np1) - T(:,n))/(theta*dtime(np1)) - &
          (one - theta)*Tdot(:,n)/theta
     ! fill energy content arrays for next time step
     if(abs(time - time_fin)/dtime(np1) < small) exit
  end do     ! end of time step loop

end program abc
