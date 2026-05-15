subroutine AO_to_MO(nBas,c,ERI_AO,ERI_MO )

! Expression of bi-electronic integrals in the MO basis set

  implicit none
!  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  double precision,allocatable  :: scr(:,:,:,:)
  double precision              :: temp1(nBas,nBas,nBas,nBas), temp2(nBas,nBas,nBas,nBas),temp3(nBas,nBas,nBas,nBas)
  !double precision,allocatable  :: scr_K(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)
  !double precision,intent(out)  :: ERI_MO_K(nBas,nBas,nBas,nBas)

! Memory allocation

allocate (scr(nBas,nBas,nBas,nBas))
!allocate (scr_K(nBas,nBas,nBas,nBas))



!--------------------------------------
! AO to MO transformation starts here !
!--------------------------------------

 write(*,*)"This is AO_to_MO"
!  write(*,*) "nBas",nBas
!write(*,*) mu,nu,la,si
!"c", c
!write(*,*)"ERI AO",ERI_AO(nBas,nBas,nBas,nBas)

!scr(nBas,nBas,nBas,nBas) =  ERI_AO(nBas,nBas,nBas,nBas)
!scr_K(nBas,nBas,nBas,nBas) =  ERI_AO(nBas,nBas,nBas,nBas)

ERI_MO = 0d0   ! set all elements to 0
!ERI_MO_K = 0d0   ! set all elements to 0 
!write(*,*)"ERI:", ERI_AO(1,1,1,:) 
 ! write(*,*)ERI_MO_J
 !  write(*,*)"ERI MO_K:"
 ! write(*,*)ERI_MO_K
 ! write(*,*)""

!do p = 1, nBas
!  do q = 1, nBas
!    do r = 1, nBas
!      do s = 1, nBas
!        do mu = 1, nBas
!          do nu = 1, nBas
!            do la = 1, nBas
!              do si = 1, nBas
              !ERI_MO(p,r,q,s)  = ERI_MO(p,r,q,s) + ( c(mu,p) * c(nu,q) * c(la,r) *c(si,s) * scr(mu, la, nu, si) )  
!               ERI_MO(p,r,q,s)  = ERI_MO(p,r,q,s) + ( c(mu,p) * c(nu,q) * c(la,r) *c(si,s) * ERI_AO(mu, la, nu, si) )
              !ERI_MO(p,q,r,s)  = ERI_MO(p,q,r,s) + ( c(p,mu) * c(q,nu) * c(r,la) *c(s,si) * scr(mu, nu, la, si) )
                !ERI_MO_K(p,q,s,r)  = ERI_MO_K(p,q,s,r) + ( c(mu,p) * c(nu,q) * c(la,s) *c(si,r) * scr_K(mu, nu, si, la) )
      !write(*,*)"c(mu,p)",c(mu,p)
      !write(*,*)"c(nu,q)",c(nu,q)
      !write(*,*)"c(la,r)",c(la,r)
      !write(*,*)"c(si,s)",c(si,s)
      !write(*,*)"scr(mu, nu, la, si)",scr(mu, nu, la, si)
!               enddo
!            enddo
!          enddo
!         enddo
!       enddo
!    enddo
!  enddo
!enddo
temp1 = 0d0
temp2 = 0d0
temp3 = 0d0

 !ERI_MO(p,r,q,s)  = ERI_MO(p,r,q,s) + ( c(mu,p) * c(nu,q) * c(la,r) *c(si,s) * ERI_AO(mu, la, nu, si) )

do p = 1, nBas
    do mu = 1, nBas

      temp1(p,:,:,:) = temp1(p,:,:,:) +   c(mu,p) * ERI_AO(mu, :, :, :) 
    enddo

    do r = 1, nBas
       do la = 1, nBas
        temp2(p,r,:,:) =  temp2(p,r,:,:) +  c(la,r) * temp1(p,la , :, :)
       enddo

       do q = 1, nBas
         do nu = 1, nBas
          temp3(p,r,q,:) =  temp3(p,r,q,:) + c(nu,q) *  temp2(p,r , nu, :)  
         enddo
        
       do s = 1, nBas
         do si = 1, nBas
           ERI_MO(p,r,q,s) = ERI_MO(p,r,q,s) + c(si,s) * temp3(p, r , q, si)
         enddo 
 enddo
 enddo
 enddo
 enddo

!write(*,*)"ERI MO_J:"
!write(*,*)"ERI_MO",ERI_MO(1,1,1,:)
!write(*,*)"ERI_MO",ERI_MO
!write(*,*)"ERI MO_K:"
! write(*,*)ERI_MO_K
! write(*,*)""

write(*,*)"This is AO_to_MO End"






end subroutine AO_to_MO
