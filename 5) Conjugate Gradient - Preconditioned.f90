program CG
   implicit none
   integer :: i,j,k,n,giri
   Real*8 :: alfa, beta, conv, tau, gamma, gammaold
   Real*8, dimension(:,:), allocatable :: a
   Real*8, dimension(:), allocatable :: x,x1,b,r,r1,h,d,z,z1
!
   write(*,*)'Inserisci dimensione sistema: '
   read(*,*)n
   allocate (a(n,n),x(n),x1(n),b(n),r(n),r1(n),h(n),d(n),z(n),z1(n))
!
!  preamble
!
   call Pream (i,j,n,conv,a,x,b,r,giri,d,beta)
!
!  Defining initial residual
!
   call dgemv('N',n,n,1.d0,a,n,x,1,0.d0,r1,1)
   r1=r1-b
!
!  Preconditioning
!
   do i=1,n
      z1(i)=r1(i)/d(i)
   end do
   gammaold = dot_product(r1,z1)
!
   do k=1,n
      giri=giri+1
!
!     Alpha term creation
!
      call dgemv('N',n,n,1.d0,a,n,z1,1,0.d0,h,1)
      gamma=dot_product(z1,h)
      alfa=-(gammaold/gamma)
!
!     creazione termine x1
!
      x1=x + alfa*(z1+beta*z)
      z=z1                       !saving z1 term
!
!     Calculating r1 and z1
!
      call dgemv('N',n,n,1.d0,a,n,x1,1,0.d0,h,1)
      r1=h-b
!
!     Preconditioning
!
      do i=1,n
         z1(i)=r1(i)/d(i)
      end do
!
!     calculating tau residual
!
      tau=norm2(r1)
!
      if(tau .lt. conv)then
         exit
      end if
!
!     Beta term creation
!
      gamma=dot_product(r1,z1)
      beta=gamma/gammaold
!
!     saving variables
!
      gammaold=gamma
      x=x1
   end do
!
   write(*,*)'Number of iterations: ', giri
   write(*,*)x
!
end program
!
!preamble
!
subroutine Pream (i,j,n,conv,a,x,b,r,giri,d,beta)
   implicit none 
   integer :: i,j,n,giri
   Real*8 :: conv,beta
   Real*8 :: a(n,n),x(n),b(n),r(n),d(n)
!
!  convergence factor
!
   write(*,*)'Insert convergence factor: '
   read(*,*)conv
!
!  Creating a matrix
!
   call SuperMatrixMaker (i,j,n,a,d)
!
!  creating vector b (known terms)
!
   do i=1,n
      b(i)=1/real(i,8) + 3.d0 
   end do
!
!  creating x vector
!
   x(i)=0.d0
!
!  Variable initialization
!
   r=0.d0
   giri=0
   beta=0.d0
end subroutine
!
!creatring Matrix a
!
subroutine SuperMatrixMaker (i,j,n,a,d)
   implicit none
   Integer :: i, j, n
   Real*8 :: a(n,n), d(n)
   do i=1,n
      a(i,i)=i+1.d0
      d(i)=a(i,i)
      do j=1,i-1                      
         a(j,i)=1.d0/(real(j,8)+real(i,8))
         a(i,j)=a(j,i)
      end do
   end do  
end subroutine
