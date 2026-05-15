subroutine CIS(nBas,nO, ERI_MO, e )

        implicit none
        include 'parameters.h'

        ! INput variables
         integer,intent(in)           :: nBas
         integer,intent(in)           :: nO
         double precision,intent(in)  :: ERI_MO(nBas,nBas,nBas,nBas)
         !double precision,intent(in)  :: ERI_MO_K(nBas,nBas,nBas,nBas)
         double precision,intent(in)  :: e(nBas)


        ! Local variables
        double precision,external     :: Kronecker_delta ! in utils
        integer                       :: a, b, i, j, OV, nV, index_ia, index_jb, id
        double precision              :: Gap 
        double precision              :: temp1
        double precision              :: temp2


        ! Output variables
         double precision,allocatable  :: A_matrix(:,:) ! matrix A
         double precision,allocatable  :: orbital_energy(:,:)
         double precision,allocatable  :: eigen(:)




  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          CIS calculation           |'
  write(*,*)'**************************************'
  write(*,*)

!allocate (ERI_MO_J(nBas,nBas,nBas,nBas))
!allocate (ERI_MO_K(nBas,nBas,nBas,nBas))
!allocate (e(:))
!allocate(orbital_energy(OV,OV))

nV = nBas - nO
OV = nO * nV  ! dimension of A
write(*,*)"OV is:", OV
write(*,*)"nV is:", nV
write(*,*)"nO is:", nO

allocate(orbital_energy(OV,OV), A_matrix(OV,OV), eigen(OV))

!   try see HOMO-LUMO gap
!    write(*,*)"e is:", e
!    write(*,*)"nBas is:", nBas
!      write(*,*)"e(nO+1):",e(nO+1) ! first virtual is nO +1
!      write(*,*)"e(nO-1):",e(nO-1)
!      write(*,*)
!      write(*,*)"e(nO):",e(nO)
!      write(*,*)"e(nV):",e(nV)
!      write(*,*)
!      write(*,*)"e(nV+1):",e(nV+1)
!      write(*,*)"e(nV-1):",e(nV-1)
!      Gap = e(nO+1) - e(nO)
!       write(*,*)
!      write(*,*)"Gap is: e(nO+1) - e(nO)",Gap
!      Gap = e(nV) - e(nO)
!       write(*,*)"Gap is:e(nV) - e(nO)",Gap 
! ****************** !

! write(*,*)"before A"
  A_matrix = 0d0
 ! write(*,*)"after A"
 !do id = 1, OV
do i = 1, nO    ! number of occupied
    do j = 1, nO
      do a = nO+1, nBas
        do b = nO+1, nBas  ! number of virtual
           if (i==j .AND. a==b) then
         !    write(*,*)"before iajb complete"
             index_ia = i + (nO * (a-nO-1)) !i * nV+ a - nBas
             !index_ia = id
!              write(*,*)"index_ia:",index_ia
             index_jb = j + (nO * (b-nO-1)) !j * nV+ a - nBas ! (j * (nV-1)) + b -nBas +1
             !index_jb = id
!             write(*,*)"index_jb:",index_jb
             !orbital_energy(index_ia, index_jb) = e(a) - e(i) ! a is from 1 to nV, i--1 to nO
              orbital_energy(index_ia, index_jb) = e(a) - e(i)
             !write(*,*)"orbital_energy(index_ia, index_jb):",orbital_energy(index_ia, index_jb)
              write(*,*)""
           !  write(*,*)"orbital energy complete"
             end if
          enddo
        enddo
      enddo
    enddo
  !  enddo

!    write(*,*)""
!write(*,*)"orbital energy:",orbital_energy
!write(*,*)""

!do j = 1,OV
!            write(*,'(1000f16.8)') orbital_energy(j,:)
!      end do

  !write(*,*)"ERI_MO_J(i,j,a,b)",ERI_MO_J
  !write(*,*)"ERI_MO_K(i,j,b,a)",ERI_MO_K
 ! write(*,*)"2* ERI_MO_J(i,j,a,b)",2* ERI_MO_J
  do i = 1, nO    ! number of occupied
    do j = 1, nO
      do a = nO+1, nBas
        do b = nO+1, nBas  ! number of virtual
         index_ia =  i + (nO * (a-nO-1)) ! i * nV+ a - nBas
         !write(*,*)"index_ia is:", index_ia
         index_jb =   j + (nO * (b-nO-1)) ! j * nV+ a - nBas
         !write(*,*)"index_jb is:", index_jb
           if (i==j .AND. a==b) then
            
             !index_ia = (i *( nV-1)) + a
             !index_ia = (i *( nV-1)) + a -nBas +1
             !index_ia = i * nV+ a - nBas
            !  write(*,*)"index_ia:",index_ia
             !index_jb = (j * (nV-1)) + b
             !index_jb = (j * (nV-1)) + b -nBas +1
              !index_jb = i * nV+ a - nBas
           !   write(*,*)"index_jb:",index_jb
            !A_matrix(index_ia, index_jb) = orbital_energy(index_ia, index_jb) + 2d0 * ERI_MO(i,j,a,b) - ERI_MO(i,j,b,a)
            A_matrix(index_ia, index_jb) = orbital_energy(index_ia, index_jb) + 2d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)
    else
            A_matrix(index_ia, index_jb) =  2d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)
            !write(*,*)"A is:", A_matrix(index_ia, index_jb)

                   !temp2 = 2 * ERI_MO_J(i,j,a,b) - ERI_MO_K(i,j,b,a)
        
                   !write(*,*) "temp2:"
        !temp1 = e *  Kronecker_delta(i,j) *  Kronecker_delta(a,b)
        !write(*,*) e *  Kronecker_delta(i,j) *  Kronecker_delta(a,b)
  !A_matrix(OV, OV) = A_matrix(OV, OV) + temp1 + temp2! (e(i) *  Kronecker_delta(i,j) *  Kronecker_delta(a,b) ) + !2 * ERI_MO_J(i,j,a,b) - ERI_MO_K(i,j,b,a) 
            end if
          enddo
        enddo
      enddo
    enddo

    ! then diagonalize A? because if Kdelta = 0 still have value
!call diagonalize_matrix(nBas, cp, e) !E is vector
  !  diagonalize (A)
!  write(*,*)""
!  write(*,*)"A matrix is:"
  !write(*,*) A_matrix
  !write(*,*) A
!  do j = 1,OV
!            write(*,'(1000f16.8)') A_matrix(j,:)
!      end do

call diagonalize_matrix(OV, A_matrix, eigen) !E is vector

!write(*,*)""
!  write(*,*)"A matrix after diagonalize is:"
!   do j = 1,OV
            !write(*,'(1000f16.8)') A_matrix(j,:)
!      end do

  !write(*,*)"The value of CIS are:"
  !    do j = 1,OV
  !          write(*,*) eigen(j)
  !    end do
  !write(*,*) eigen

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'              CIS (excitation) energy                            '
  write(*,*)'----------------------------------------------------'
  write(*,*)' E(CIS) = ' 
  !write(*,*) !eigen(1)
             do j = 1,OV
            write(*,*) eigen(j)
      end do
  write(*,*)'----------------------------------------------------'
  write(*,*)

end subroutine CIS

