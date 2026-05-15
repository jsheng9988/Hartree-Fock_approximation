subroutine TDHF(nBas,nO, ERI_MO, e )

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
         double precision,allocatable  :: A_matrix(:,:) , B_matrix(:,:), AB(:,:), temp_AB(:,:)! matrix A
         double precision,allocatable  :: orbital_energy(:,:)
         double precision,allocatable  :: eigen(:)
         double precision,allocatable  :: C(:,:), X(:,:), eAB(:), e_C(:)




  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          TDHF calculation           |'
  write(*,*)'**************************************'
  write(*,*)


nV = nBas - nO
OV = nO * nV  ! dimension of A
write(*,*)"OV is:", OV
write(*,*)"nV is:", nV
write(*,*)"nO is:", nO

allocate(orbital_energy(OV,OV), A_matrix(OV,OV), B_matrix(OV, OV), eigen(OV), AB(OV,OV), temp_AB(OV,OV))
allocate(C(OV, OV), X(OV, OV), eAB(OV), e_C(OV))

! ****************** !

  A_matrix = 0d0
do i = 1, nO    ! number of occupied
    do j = 1, nO
      do a = nO+1, nBas
        do b = nO+1, nBas  ! number of virtual
           if (i==j .AND. a==b) then
         !    write(*,*)"before iajb complete"
             index_ia = i + (nO * (a-nO-1))
             !index_ia = id
             ! write(*,*)"index_ia:",index_ia
             index_jb = j + (nO * (b-nO-1)) ! (j * (nV-1)) + b -nBas +1
             !index_jb = id
             !write(*,*)"index_jb:",index_jb
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

do j = 1,OV
!            write(*,'(1000f16.8)') orbital_energy(j,:)
      end do

  !write(*,*)"ERI_MO_K(i,j,b,a)",ERI_MO_K
 ! write(*,*)"2* ERI_MO_J(i,j,a,b)",2* ERI_MO_J
  do i = 1, nO    ! number of occupied
    do j = 1, nO
      do a = nO+1, nBas
        do b = nO+1, nBas  ! number of virtual
           !if (i==j .AND. a==b) then
            
             index_ia = i + (nO * (a-nO-1)) !i * nV+ a - nBas
            !  write(*,*)"index_ia:",index_ia
              index_jb =  j + (nO * (b-nO-1))
           !   write(*,*)"index_jb:",index_jb
            !A_matrix(index_ia, index_jb) = orbital_energy(index_ia, index_jb) + 2d0 * ERI_MO(i,j,a,b) - ERI_MO(i,j,b,a)
            !A_matrix(index_ia, index_jb) = orbital_energy(index_ia, index_jb) + 2d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)
            if (i==j .AND. a==b) then
            A_matrix(index_ia, index_jb) = orbital_energy(index_ia, index_jb) + 2d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)
    else
            A_matrix(index_ia, index_jb) =  2d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)

            !write(*,*)"A is:", A_matrix(index_ia, index_jb)

        
            end if
          enddo
        enddo
      enddo
    enddo

    ! then diagonalize A? because if Kdelta = 0 still have value
!call diagonalize_matrix(nBas, cp, e) !E is vector
  !  diagonalize (A)
 ! write(*,*)""
!  write(*,*)"A matrix is:"
  !write(*,*) A_matrix
  !write(*,*) A
  do j = 1,OV
!            write(*,'(1000f16.8)') A_matrix(j,:)
      end do

!call diagonalize_matrix(OV, A_matrix, eigen) !E is vector

!write(*,*)""
!  write(*,*)"A matrix after diagonalize is:"
   do j = 1,OV
           !write(*,'(1000f16.8)') A_matrix(j,:)
      end do

!write(*,*)"eigenvalue is:"
!  write(*,*) eigen

B_matrix = 0d0
! B_matrix
 do i = 1, nO    ! number of occupied
    do j = 1, nO
      do a = nO+1, nBas
        do b = nO+1, nBas  ! number of virtual
           if (i==j .AND. a==b) then

             index_ia = i + (nO * (a-nO-1))
            !  write(*,*)"index_ia:",index_ia
              index_jb = j + (nO * (b-nO-1))
           !   write(*,*)"index_jb:",index_jb
           ! B_matrix(index_ia, index_jb) = 2 * ERI_MO(i,j,a,b) - ERI_MO(i,j,b,a)
           B_matrix(index_ia, index_jb) = 2 * ERI_MO(i,j,a,b) - ERI_MO(i,j,b,a) 
!           write(*,*)"B is:", B_matrix(index_ia, index_jb)


            end if
          enddo
        enddo
      enddo
    enddo


!     write(*,*)""
!  write(*,*)"B matrix is:"
  !write(*,*) A_matrix
  !write(*,*) A
!  do j = 1,OV
!            write(*,'(1000f16.8)') B_matrix(j,:)
!      end do



! calculate the (A - B) ^ (1/2)
temp_AB = A_matrix - B_matrix
!temp_AB = sqrt(temp_AB)
!write(*,*)"squarooted"
!write(*,*)"A-B matrix is:"
!do j = 1,OV
!            write(*,'(1000f16.8)') temp_AB(j,:)
!     end do

call diagonalize_matrix(OV,temp_AB,eAB) 


!write(*,*)"eAB before square",eAB

do i=1,OV
  eAB(i) = sqrt(eAB(i))  ! --calculate diagonal ^ 1/2
enddo

!write(*,*)""
!write(*,*)"eAB after square",eAB

call ADAt(OV, temp_AB , eAB , X)

AB = A_matrix + B_matrix
!write(*,*)"A + B matrix is:"
!do j = 1,OV
!            write(*,'(1000f16.8)') AB(j,:)
!      end do

!C = matmul( X ,matmul( AB, X))
C =  X * AB * X
!C = temp_AB * AB * temp_AB

!write(*,*)"C before diagonalize:"
!do j = 1,OV
!            write(*,'(1000f16.8)') C(j,:)
!      end do

! diagonalize C
call diagonalize_matrix(OV,C,e_C)

!write(*,*)"C after diagonalize:"
!do j = 1,OV
!            write(*,'(1000f16.8)') C(j,:)
!      end do


!write(*,*)""
!write(*,*)"e_C is",e_C

!e_C = e_C**(0.5)
!write(*,*)"Values of TDHF are:"
!do j = 1,OV
!       write(*,*) e_c(j)
!end do
!write(*,*)"Values of TDHF are:"
!do j = 1,OV
!       write(*,*) e_C(j)
!end do

write(*,*)

e_C = sqrt(e_C)
!write(*,*)"Values of TDHF are:"
!do j = 1,OV
!       write(*,*) e_C(j)
!end do

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'              TDHF energy                            '
  write(*,*)'----------------------------------------------------'
  write(*,*)' E(TDHF) = ' 
  !write(*,*) !e_C(1)
  do j = 1,OV
       write(*,*) e_C(j)
end do
  write(*,*)'----------------------------------------------------'
  write(*,*)

end subroutine TDHF

