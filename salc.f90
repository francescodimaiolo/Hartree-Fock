program salc
  implicit none
  integer:: i,b,c,dim,a, binarysearch,nso,nsiti,temp,temp2,j,m
  integer,allocatable:: uso(:), func(:)
  logical::bool, bool1
  external binarysearch

  open(1,file='configurazioni.dat')
  open(2,file='dim.dat')
  open(3,file='analisi.dat')
  open(4,file='salc.dat')
  read(2,*) dim
  allocate(uso(dim), func(dim))
  nso=26
  nsiti=13
  do i=1,dim
     read(1,*) a
     func(i)=a
  enddo
  uso=1
  do i=1,dim
     if(uso(i)==0)goto 300
     do j=0,12,2
        bool=btest(func(i),j)
        bool1=btest(func(i),j+1)
        if(bool)then
           b=1
        else
           b=0
        endif

        if(bool1)then
           c=2
        else
           c=0
        endif

        if((b+c).eq.1)then
           temp=ibclr(func(i),j)
           temp=ibset(temp,j+1)
        endif

        if((b+c).eq.2)then
           temp=ibclr(func(i),j+1)
           temp=ibset(temp,j)
        endif
     enddo

     do j=nsiti+1,nso-2,2
        bool=btest(func(i),j)
        bool1=btest(func(i),j+1)
        if(bool)then
           b=1
        else
           b=0
        endif

        if(bool1)then
           c=2
        else
           c=0
        endif

        if((b+c).eq.1)then
           temp2=ibclr(temp,j)
           temp2=ibset(temp2,j+1)
        endif

        if((b+c).eq.2)then
           temp2=ibclr(temp,j+1)
           temp2=ibset(temp2,j)
        endif
        m=binarysearch(1,dim,func,temp2)
        if(m.gt.i)then
           uso(i)=1
           uso(m)=0
        endif
     enddo
300  continue
  enddo

  do i=1,dim
     write(3,*) func(i), uso(i)
     if(uso(i).eq.1)then
        write(4,*) func(i)
     endif
  enddo
  
end program salc
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
