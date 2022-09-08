program hamiltonianci

  implicit none

  integer,allocatable:: config(:)
  integer::i,ii,jj,kk,ll, dim,n,nso, j,k, a, c, temp, m, binarysearch, nsiti, p,q,&
       x, l,n1,lda,lwork,info,junk,count1,nn,conta,nexc,alpha,beta,gamma,delta

  real*8:: energy, b, t, pot, fock, uc, un, esc,esnpy,esnaza,sz,spin
  real*8, allocatable:: ham(:,:) , emo(:), v(:,:), w(:), work(:), f(:,:),&
       xvec(:), yvec(:), dist(:,:), ea(:), ai(:), uvec(:), esite(:), vmat(:,:,:,:)

  logical:: bool,bool1,bool2,bool3,booli,boolj,boolk,booll
  character*1::jobz,uplo
 
  external binarysearch

  open(1,file='configurazioni.dat')
  open(2,file='dim.dat')
  open(3,file='autovettori-hf.dat')
  open(4,file='hamiltonianci.dat')
  open(5,file='autovalori-hf.dat')
  open(6,file='interatomic-distances.dat')
  open(7,file='input.dat')
  open(9,file='autovalori_ci.dat')
  open(8,file='spin.dat')
  open(10,file='nexc.dat')

  nsiti=13
  nso=2*nsiti
  
  Uc=11.26d0
  Un=15d0
  esc=0d0
  esnpy=-12d0
  esnaza=-2d0
  t=1

  read(2,*) dim
  close(2)

  read(8,*)spin
  close(8)

  read(10,*)nexc
  close(10)

  allocate(config(dim), ham(dim,dim), emo(nsiti),v(nsiti,nsiti),&
       f(nsiti,nsiti),yvec(nsiti),xvec(nsiti),dist(nsiti,nsiti),&
       ea(nsiti), ai(nsiti), uvec(nsiti), esite(nsiti),vmat(nsiti,nsiti,nsiti,nsiti))

  !!! reading basis functions
  do i=1,dim
     read(1,*) a
     config(i)=a
  enddo
  close(1)

  !!! reading HF MOs energies (i.e., emo)
  do i=1,nsiti
     read(5,*)junk,emo(i)
  enddo
  close(5)

  !!! reading geometry
  do i=1,13
     read(6,173)(dist(i,j),j=1,nsiti)
  enddo
  close(6)

173 format(<nsiti>(2x,f10.5))

  !! making auxiliary arrays containing U and site energies
  do i=1,nsiti
     uvec(i)=Uc
     esite(i)=esc
     if(i.eq.nsiti)then
        uvec(i)=Un
        esite(i)=esnpy
     endif
  enddo

  !! getting the electrostatic potential
  do i=1,nsiti-1
     v(i,i)=uvec(i)
     do p=i+1,nsiti
        v(i,p) = 14.397d0 &
             / dsqrt( dist(i,p)**2 + (28.794d0/(uvec(i)+uvec(p)))**2 )
        V(p,i)=V(i,p)
     enddo
  enddo
  v(nsiti,nsiti)=uvec(nsiti)

!!! reading HF MOs
  do i=1,nsiti
     do j=1,nsiti
        read(3,*) a,c,fock
        f(a,c)=fock
     enddo
  enddo
  close(3)


  vmat=0d0
  do alpha=1,nsiti
     do beta=1,nsiti
        do gamma=1,nsiti
           do delta=1,nsiti
              
              do ii=1,nsiti
                 do jj=1,nsiti
                    vmat(alpha,beta,gamma,delta) = vmat(alpha,beta,gamma,delta) + &
                         f(ii,alpha)*f(ii,beta)*f(jj,gamma)*f(jj,delta)*v(ii,jj)
                 enddo
              enddo

           enddo
        enddo
     enddo
  enddo
  
  !=========================DIAGONALE=========================
  
  ham=0d0
  do n=1,dim
     energy=0.d0
     do i=0,nso-2,2
        bool=btest(config(n),i)
        bool1=btest(config(n),i+1)
        if(bool)then
           energy=energy+emo((i+2)/2)
        endif
        if(bool1)then
           energy=energy+emo((i+2)/2)
        endif
     enddo
     ham(n,n)= energy
     write(89,*) n,ham(n,n)
  enddo

  !===========================================================
  !=========================FUORI DIAGONALE=========================
  do n=2,dim  !! we do not account for the GS (i.e., n==1)

     !if(n.eq.3)write(*,*)i,j,k
     
     do i=0,nso-1
        alpha=i/2+1
        do j=0,nso-1
           beta=j/2+1
           do k=0,nso-1
              gamma=k/2+1
              do l=0,nso-1
                 delta=l/2+1

                 bool3=btest(config(n),k)
                 if(bool3)then
                    count1=0
                    temp=ibclr(config(n),k)
                    do kk=0,k
                       boolk=btest(temp,kk)
                       if(boolk)count1=count1+1
                    enddo

                    bool2=btest(temp,l)
                    if(bool2)then
                       temp=ibclr(temp,l)
                       do ll=0,l
                          booll=btest(temp,ll)
                          if(booll)count1=count1+1
                       enddo

                       bool1=btest(temp,j)
                       if((.not.bool1) .and. &
                            ( ( ((j/2)*2.eq.j) .and. ((l/2)*2.eq.l) ) .or. &
                              ( ((j/2)*2.ne.j) .and. ((l/2)*2.ne.l) ) ) )then
                          do jj=0,j
                             boolj=btest(temp,jj)
                             if(boolj)count1=count1+1
                          enddo
                          temp=ibset(temp,j)

                          bool=btest(temp,k)
                          if((.not.bool) .and. &
                               ( ( ((i/2)*2.eq.i) .and. ((k/2)*2.eq.k) ) .or. &
                                 ( ((i/2)*2.ne.i) .and. ((k/2)*2.ne.k) ) ) )then
                             do ii=0,i
                                booli=btest(temp,ii)
                                if(booli)count1=count1+1
                             enddo
                             temp=ibset(temp,i)
                             
                             sz=0d0
                             do ii=0,nso-2,2
                                bool=btest(temp,ii)
                                bool1=btest(temp,ii+1)
                                if(bool)sz=sz+0.5d0
                                if(bool1)sz=sz-0.5d0
                             enddo

                             conta=0
                             do ii=nsiti+1,nso-1
                                bool=btest(temp,ii)
                                if(bool)conta = conta + 1
                             enddo

                             if((dabs(sz-spin).lt.1d-8) .and. (conta.eq.nexc))then
                                m = binarysearch(1,dim,config,temp)
                                if(n.eq.79.and.m.eq.65)&
                                !if((n.eq.85).and.(m.eq.80))&
                                     write(87,'(6i3,f10.5)')n,m,i,j,l,k,((-1d0)**count1)
                                if(n.ne.m)ham(n,m)=ham(n,m)+&
                                     0.5d0*((-1d0)**count1)*vmat(alpha,beta,gamma,delta)
                             endif
                             
                          endif
                       endif
                    endif
                 endif
                 
              enddo
           enddo
        enddo
        
     enddo

  enddo
     
  do i=1,dim
     do j=1,dim
        if(dabs(ham(i,j)).gt.1d-8)write(4,*)i,j,ham(i,j)
        !write(4,'(<dim>(f10.5))')(ham(i,j),j=1,dim)
     enddo
  enddo
     
  jobz='V'
  uplo='U'
  n1=dim
  lda=dim
  lwork=3*n-1
  
  allocate(w(n),work(lwork))
  call dsyev (jobz,uplo,n1,ham,lda,w,work,lwork,info)
  write(*,*) info
  do i=1,dim
     write(9,*) w(i)
  enddo
end program

!=========================BINARY SEARCH=========================

integer function binarysearch(i, length, array, val)
  ! Given an array and a value, returns the index of the element that
  ! is closest to, but less than, the given value.
  ! Uses a binary search algorithm.
  ! "delta" is the tolerance used to determine if two values are equal
  ! if ( abs(x1 - x2) <= delta) then
  ! assume x1 = x2
  ! endif

  implicit none

  integer, intent(in) :: length, i
  integer, dimension(length), intent(in) :: array
  integer, intent(in) :: val

  !integer :: binarysearch

  integer :: left, middle, right

  left = i
  right = length
  binarysearch=0


  if (val.lt.array(left) .or. val.gt.array(right)) go to 10

  do

     if (left .gt. right) then
        exit
        !write(*,*) 'ERRORE!!!'
     endif

     !divisione=((left+right) / 2.0)
     !middle = jnint(divisione)
     middle=(left+right)/2

     if ( array(middle) .eq. val ) then
        binarySearch = middle
        return
     else 
        if (array(middle) .gt. val) then
           right = middle - 1
        else
           left = middle + 1
        end if
     end if
  end do

  binarysearch = right
10 continue
end function binarysearch
!==================================================

