program hartree
  implicit none

  integer::i,iii,j,k,m,l, n, p, nsiti, lda, lwork, info,niter,valence(13), aaa
  integer, allocatable :: nz(:)
  
  real*8:: n1, xvec(13), yvec(13), V(13,13), opj(13,13), Un, Un2, Uc, V1, V2, opk(13,13),&
       k0, j0, a1, b, a2, jj, opj2(13,13), opk2(13,13), g(13,13), h(13,13), ii(13,13),&
       opj3(13,13),opk3(13,13),fock3(13,13),fock4(13,13), t, energia, opj4(13,13), opk4(13,13),&
       thresh,eold, kappa, alpha, l1, pi,x,y, ll, sum,tcn,esc,esnpy,esnaza
  real*8, allocatable:: ham(:,:), a(:,:), w(:), work(:), f(:,:), huck(:,:),work1(:), w1(:),&
       work2(:), w2(:), U(:), es(:), dist(:,:), ai(:), ea(:), denmat(:,:), denmatold(:,:),&
       uvec(:)
  
  character*1::jobz, uplo,labels(13)

  
  open(1,file='huckel1.dat')
  open(2,file='autovettori-huckel1.dat')
  open(3,file='potenziali.dat')
  open(4,file='operatore-j.dat')
  open(5,file='operatore-k.dat')
  open(6,file='fock-1.dat')
  open(7,file='input.dat')
  open(8,file='interatomic-distances.dat')
  open(20,file='autovalori-huckel.dat')
  open(21,file='autovalori-fock.dat')
  open(112,file='eigenvec_HF.dat')
  open(114,file='check.dat')

  nsiti=13

  !!! Model parameters
  t=-2.4  !! C-C hopping integral
  tcn=-2.5d0  !! C-N hopping integral
  Uc=11.26d0  !! C hubbard U
  Un=15.d0  !! N hubbard U
  esc=0d0  !! C site energy (taken as zero)
  esnpy=-12d0 !! pyrrole N site energy
  esnaza=-2d0 !! aza N site energy

  !!! Allocation
  allocate(ham(nsiti,nsiti), f(nsiti,nsiti),huck(nsiti,nsiti),U(nsiti),&
       nz(nsiti),es(nsiti),dist(nsiti,nsiti),ai(nsiti),ea(nsiti),&
       denmat(nsiti,nsiti),denmatold(nsiti,nsiti),uvec(nsiti))

  !!! Creating auxiliary arrays
  do i=1,nsiti
     uvec(i)=Uc
     es(i)=esc
     nz(i)=1
     if(i.eq.13)then
        uvec(i)=Un
        es(i)=esnpy
        nz(i)=2
     endif
     write(114,'(2(i3),2(f10.5))')i,nz(i),es(i),uvec(i) !! just for checking purpose
  enddo

  close(114)

!!! Writing the starting hamiltonian !!!

  ham=0d0
  do i=1,nsiti-1
     
     ham(i,i)=es(i)
     ham(i,i+1)=t
         
     if(i==1)ham(i,12)=t !! PBC
     if(i==4)ham(i,13)=t
     if(i==8)ham(i,13)=t
        
  enddo
  ham(nsiti,nsiti)=es(nsiti)

  do i=1,nsiti
     do j=1,nsiti
        ham(j,i)=ham(i,j)
     enddo
  enddo

  do i=1,nsiti
     write(1,'(<nsiti>(f10.5))')(ham(i,j),j=1,nsiti)
  enddo

  close(1)
  
  huck=ham !! making a copy of the starting hamiltonian
  
  jobz='V'
  uplo='U'
  n=nsiti
  lda=nsiti
  lwork=3*N-1
  allocate(w(n), work(lwork))

  call dsyev (jobz,uplo,n,ham,lda,w,work,lwork,info)
  write(*,*) 'info=', info

  do i=1,nsiti
     write(20,*) w(i)
  enddo

  energia=0.d0
  do i=1,7
     energia= energia+2*w(i)
  enddo
  write(20,*) 'GS energy=',energia

  !! writing the density matrix
  denmat=0d0
  do i=1,nsiti
     do j=1,nsiti

        do k=1,7 !! loop on occupied MOs
           denmat(i,j) = denmat(i,j) + 2d0 * ham(i,k) * ham(j,k)
        enddo
        
     enddo
  enddo
  
  do i=1,nsiti
     write(2,'(<nsiti>(f10.4))')(ham(i,j),j=1,nsiti) 
  enddo

  close(2)
  
  !==================================================

  !!! molecular geometry
  l1=1.50
  pi=dacos(-1.d0)
  alpha=pi/6
  do i=1,13
     if(i==1)then
        x=-l1*dcos(alpha)
        y=l1*(1+dsin(alpha))
     endif
     if(i==2)then
        x=-2*l1*dcos(alpha)
        y=l1
     endif
     if(i==3)then
        x=-2*l1*dcos(alpha)
        y=0
     endif
     if(i==4)then
        x=-l1*dcos(alpha)
        y=-l1*dsin(alpha)
     endif
     if(i==5)then
        x=-l1*dcos(alpha)
        y=-l1*(dsin(alpha)+1)
     endif
     if(i==6)then
        x=0
        y=-l1*(2*dsin(alpha)+1)
     endif
     if(i==7)then
        x=l1*dcos(alpha)
        y=-l1*(dsin(alpha)+1)
     endif
     if(i==8)then
        x=l1*dcos(alpha)
        y=-l1*dsin(alpha)
     endif
     if(i==9)then
        x=2*l1*dcos(alpha)
        y=0
     endif
     if(i==10)then
        x=2*l1*dcos(alpha)
        y=l1
     endif
     if(i==11)then
        x=l1*dcos(alpha)
        y=l1*(1+dsin(alpha))
     endif
     if(i==12)then
        x=0
        y=l1
     endif
     if(i==13)then
        x=0
        y=0
     endif
     xvec(i)=x
     yvec(i)=y
     write(12,*) x, y, i
  enddo
  do i=1,nsiti
     do p=1,nsiti
        ll=((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2)**0.5
     !   write(*,*) ll
        dist(i,p)=ll
     enddo
  enddo
  do i=1,nsiti
     write(8,'(<nsiti>(2x,f10.5))')(dist(i,j),j=1,nsiti)
  enddo

  !==================================================
  
  do i=1,nsiti-1
     v(i,i)=uvec(i)
     do p=i+1,nsiti
        v(i,p) = 14.397d0 &
             / dsqrt( ((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2) + (28.794d0/(uvec(i)+uvec(p)))**2 )
        v(p,i)=v(i,p)
     enddo
  enddo
  v(nsiti,nsiti)=uvec(nsiti)

  do i=1,nsiti
     write(3,'(<nsiti>(f10.4))')(v(i,j),j=1,nsiti)
  enddo

!!! Now, starting the Hartree-Fock iterations
  
  niter=10000000
  thresh=1d-9
  eold=energia
  denmatold=denmat
  
  do iii=1,niter

!!!OPERATORE j
     opj=0.d0; opk=0.d0
     do m=1,nsiti
        do l=1,nsiti
           opj(m,m) = opj(m,m) + (denmat(l,l)-nz(l)*1d0) * V(m,l)
        enddo
     enddo
     
     do i=1,nsiti
        write(4,'(<nsiti>(2x,f10.5))')(opj(i,j),j=1,nsiti)
     enddo

!!! OPERATORE K

     do m=1,nsiti
        opk(m,m) = (0.5d0*denmat(m,m)-nz(m)*1d0) * v(m,m)
        do n=1,nsiti
           if(m.ne.n)opk(m,n) = 0.5d0*denmat(n,m) * v(n,m)
        enddo
     enddo

     do i=1,nsiti
        write(5,'(<nsiti>(2x,f10.5))')(opk(i,j),j=1,nsiti)
     enddo
     !=======OPERATORE DI FOCK=========================
     f=0d0
     do i=1,nsiti
        do j=1,nsiti
           f(i,j) = huck(i,j) + opj(i,j) - opk(i,j)
        enddo
     enddo

     !!==========Diagonalizzazione fock==========

     deallocate(w,work)
     jobz='V'
     uplo='U'
     n=nsiti
     lda=nsiti
     lwork=3*n-1

     allocate(w(n),work(lwork))
     call dsyev (jobz,uplo,n,f,lda,w,work,lwork,info)

     energia=0.d0
     do i=1,7
        energia= energia+2*w(i)
     enddo
     write(21,*) 'energia=',energia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Updating the coefficient matrix

     ham=0d0
     do i=1,nsiti
        do j=1,nsiti
           ham(i,j)=f(i,j)
        enddo
     enddo

     denmat=0d0
     do i=1,nsiti
        do j=1,nsiti

           do k=1,7
              denmat(i,j) = denmat(i,j) + 2d0 * ham(i,k) * ham(j,k)
           enddo
           
        enddo
     enddo

!!! Checking SCF convergence
     do i=1,nsiti
        do j=1,nsiti
           !if(dabs(eold-energia).gt.thresh)then
           if(dabs(denmatold(i,j)-denmat(i,j)).gt.thresh)then
              eold=energia
              denmatold=denmat
              go to 129  !! keep on iterating
           endif
        enddo
     enddo
     go to 130  !! converged !!

129  continue
  enddo

130 continue

  write(*,131)'SCF convergence reached after',iii,'iterations, with threshold=',thresh

131 format (1x,a29,2x,i5,2x,a27,2x,e12.5)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

labels=''
labels(7)='*'


write(112,'(<nsiti>(f10.4))')(w(i),i=1,nsiti)
write(112,'(<nsiti>(2x,a10))')(labels(i),i=1,nsiti)
write(112,*)
do i=1,nsiti
   write(112,'(<nsiti>(f10.4))')(ham(i,j),j=1,nsiti)
   write(113,*)i,w(i)
enddo

write(*,*)'HOMO energy=',w(7),'eV'
write(*,*)'LUMO energy=',w(8),'eV'
write(*,*)'HOMO-LUMO gap=',w(8)-w(7),'eV'


end program hartree


