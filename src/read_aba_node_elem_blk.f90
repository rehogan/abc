
subroutine read_aba_node_elem_blk
  !     things related to element blocks are stored by material ID
  !     instead of the sequence in which they are read in
  !     this is done so that things are grouped by material ID
  !     instead of element block ID

  ! glossary:
  ! no_el_blk       = no of element blocks
  ! no_el_in_blk(i) = array of no of elements in block i
  ! mat_eb(i)       = material ID for element block i

  use mod_interface_definitions
  use mod_precision
  use mod_constants
  use mod_global_parameters
  use mod_nodes_elements
  use mod_files
  use mod_parsing
  implicit none
  integer                 :: nte, ne, mat_id_aba, eof, lnn
  integer                 :: i, i12, i13, i23, i24, i34, i41, j, jj, k, kk
  real(kind=dp)           :: tinit_uni

  integer, dimension(3):: nndx

  itinit = 'no'
  ncol = 60              ! when parsing "aba" file, only look at 1st 60 col
  read(kaba,"(/a72)") grid_title
  read(kaba,"(//40x,i6)") no_nodes ! skip 2 lines, read # nodes
  read(kaba,"(40x,i6)") no_elem    ! read # elements
!  read(kaba,"(40x,i6)") nte   !               read # elements
  read(kaba,"(///)")          ! skip 3 lines
  ! allocate memory based on no of total nodes (no_nodes)
  allocate( x(no_nodes), y(no_nodes), xi(no_nodes), yi(no_nodes) )
  allocate( T(no_nodes,3), Tdot(no_nodes,3), dt(no_nodes), t_li(no_nodes) )
  ! allocate energy content en_cnt(i,j,k) 
  !              i = # elem, j = # time indcees, k = # nodes/elem
  allocate( en_cnt(no_elem,3,4) )
  allocate( nn(4,no_elem))
  !     ready to read nodal data, ignore z-coordinate
  do i=1,no_nodes               ! loop over total # nodes
     read(kaba,*) j,xi(i),yi(i),zi_dum ! i in xi and yi signifies input units
     if(itinit == 'yes') then
        t(i,1) = zi_dum    ! initial temp from .rst file
        t(i,2) = zi_dum    ! initial temp from .rst file
     else
        t(i,1) = tinit_uni ! uniform initial temp from .inp file
        t(i,2) = tinit_uni ! uniform initial temp from .inp file
     end if
     x(i) = xi(i)/len_conv
     y(i) = yi(i)/len_conv
  end do
  Tdot(1:no_nodes,1:3) = zero
  !     read the element blocks

  nhb = 0
  no_el_blk = 0            ! initialize counter for no element blocks
  jj = 0                   ! sequential element number counter init
  do k=1,max_no_el_blk                ! maximum no element blocks
     read(kaba,"(a80)",iostat=eof) record  ! read element header record or eof
     if(eof == -1 .or. record(1:5) == '*NSET' .or. record(1:6) == '*ELSET') then
     
   exit        ! found end of file or end of useful data in .aba
        ! looking for "MAT" to indicate beginning of element block
     else if(compare_string(record,' = MAT')) then ! found element block header
        ncol = index(record,'MAT')
        read(record(ncol+3:ncol+5),"(i3)") mat_id_aba
        no_el_blk = no_el_blk + 1
        no_el_in_blk(k) = 0  ! init el in blk counter
        kk = 0                            ! initialize counter
     endif
     ! read connectivity for this element block
     do j=1,max_elem_in_blk
        read(kaba,"(a80)",iostat=eof) record  ! read element record
        if(eof == -1) then
           exit
        end if
        call record_type(ncol,ncol_max,record,rec_typ)
        if(rec_typ == 'c') then
           backspace(kaba)
           exit
        end if
        read(record,*) ne,(nn(lnn,ne),lnn=1,4)
        jj = jj + 1
        seq_el_no(jj) = ne
        no_el_in_blk(k) = no_el_in_blk(k) + 1
     end do
     mat_no(k) = mat_id_aba
     kk = kk + 1
     ! sequentially number elements in each element blk
     ne_eb(kk,mat_id_aba) = ne       ! list elements by blk
     if (nn(4,ne) == 0) nn(4,ne) = nn(3,ne)     ! 3 node tri
     if (mat_no(ne) == 0) then
        mat_no(ne) = mat_no(ne - 1) ! 3 node tri
     end if
     !  determine half bandwidth of matrix
     i12 = iabs(nn(1,ne) - nn(2,ne))
     i23 = iabs(nn(2,ne) - nn(3,ne))
     i34 = iabs(nn(3,ne) - nn(4,ne))
     i41 = iabs(nn(4,ne) - nn(1,ne))
     i13 = iabs(nn(1,ne) - nn(3,ne))
     i24 = iabs(nn(2,ne) - nn(4,ne))
     nhb = max0(i12,i23,i34,i41,i13,i24,nhb)
  end do


  nfb = 2*nhb + 1    ! full band width
  nhbp1 = nhb + 1    ! half band width + 1

  ! allocate arrays c, b, and z
  allocate ( c(no_nodes, 2*nhb+1), b(no_nodes+nhb+1), z(no_nodes+nhb+1) )

  ! verify that connectivity data properly read
  ! c(:,:) = global matrix in band storage format
  ! b(:)   = right hand side vector + some for subroutine band
  ! z(:)   = temporary array
  jj = 0
  do k=1,no_el_blk           ! loop on element blocks
     write(6,"('element block = ',i3,' material no = ',i3)") k, mat_no(k)
     do j=1,no_el_in_blk(k)  ! loop on # elements in element block k
        jj = jj + 1          ! sequential element no
        kk = seq_el_no(jj)
        write(6,"(6i4)") jj,kk,(nn(lnn,kk),lnn=1,4)
     end do
  end do
  write(6,"(' no_el_blk = ',i3, ' no_elem_in_blk = ', i3)") & 
       no_el_blk, no_el_in_blk(1)
  nte = sum(no_el_in_blk(1:no_el_blk))
  if(no_elem /= nte) then
     write(6,"('********check no_elem and nte*****')")
     write(6,"(' no_elem = ', i3, ' nte = ',i3)") no_elem, nte
  end if
  
end subroutine read_aba_node_elem_blk
