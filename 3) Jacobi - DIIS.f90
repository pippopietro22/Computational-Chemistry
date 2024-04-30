program DIIS
   Implicit none
   Real*8 :: errore, tau
   Real*8, dimension(:,:), allocatable :: A, O, B_mat, e_mat, x_mat 
   Real*8, dimension(:), allocatable :: x, x_new, b, diagonal, error, err_new, Coeff, T, y
   integer :: i,j,k
   Integer :: n,giri,cicli, dimB, it_max
   Integer :: info, ipiv(21)
!
!  Insert Matrix Dimension
!
   write(*,*)'Insert dimension: '
   read(*,*)n
!
   allocate(A(n,n),O(n,n), x(n), x_new(n), b(n), diagonal(n), y(n), error(n), err_new(n),&
            B_mat(21,21), x_mat(n,20), e_mat(n,20), t(21), coeff(20))
!
   call pream (i,j,n,tau,a,o,diagonal,x,b)
!
   giri=0
   Cicli=0
   it_max=200
!
!  Jacobi begins
!
   do k=1,it_max
      giri=giri+1
!  
      if(giri .eq. 20)then
         cicli=cicli+1
      end if     
!
      if(cicli .eq. 0)then
            dimB=giri
         else
            dimB=20 
      end if
!
!     Jacobi
!
      call dgemv('N',n,n,1.d0,o,n,x,1,0.d0,Y,1)
      x_new=b-y
!
      x_new=x_new/diagonal
!
!     Saving xnew
!
      x_mat(:,giri)=x_new(:)
!
!     DSSI      
!      
      err_new = x_new-x
!
!     Saving vectors err_new
!
      e_mat(:,giri)=err_new(:)
!
!     Creating matrix B
!
      B_mat=0.d0
      B_mat(:,dimB+1)=1.d0
      B_mat(dimB+1,:)=1.d0
      B_mat(dimB+1,dimB+1)=0.d0
      do i=1,dimB
         B_mat(i,i)=dot_product(e_mat(:,i),e_mat(:,i))
         do j=1,i-1
            B_mat(j,i)=dot_product(e_mat(:,i),e_mat(:,j))
            B_mat(i,j)=B_mat(j,i)
         end do
      end do
!
!     creating known terms
!
      T=0.d0
      T(dimB+1)=1.d0
!
!     Solving the "small" linear system
!
      call dgesv(dimB+1,1,B_mat,21,ipiv,T,dimB+1,info)
      if(info .ne. 0)then
         write(*,*)'dgesv è fallita con errore', info
         stop
      end if
!      
      do i=1,20
         Coeff(i)=T(i)
      end do
!
!     ricreating vectors err_new
!
      call dgemv('N',n,20,1.d0,e_mat,n,Coeff,1,0.d0,err_new,1)
!   
!     ricreating vectors x_new
!
      call dgemv('N',n,20,1.d0,x_mat,n,Coeff,1,0.d0,x_new,1)
!     
!     checking error
!
      errore= norm2(err_new)
!
      if(errore .lt. tau)then 
         write(*,*)x_new
         exit
      end if
!
      x=x_new
   end do
!
end program
!
!preamble
!
subroutine pream (i,j,n,tau,a,o,diagonal,x,b)
   implicit none
   Integer :: i,j,n
   Real*8 :: tau
   Real*8 :: a(n,n), o(n,n), diagonal(n),x(n), b(n)
!
!  Inserit Convergence Factor
!
   write(*,*)'Insert convergence Factor:'
   read(*,*)tau
!
!  Creating a matrix positive defined
!
   call SuperMatrixMaker (i,j,n,a,diagonal,o)
!
!  Creating guess vector x
!
   x=0.d0
!  Creating vector b (known term)
!
   do i=1,n
      b(i)=1/real(i,8) +3.d0
   end do
end subroutine
!
!Subroutine: Creates a matrix positive defined
!
subroutine SuperMatrixMaker (i,j,n,a,diagonal,o)
   implicit none
   Integer :: i, j, n
   Real*8 :: a(n,n),diagonal(n),o(n,n)
      do i=1,n
         a(i,i)=i+1.d0
         diagonal(i)=a(i,i)
         do j=1,i-1                      
         a(j,i)=1.d0/(real(j,8)+real(i,8))
         a(i,j)=a(j,i)
         o(i,j)=a(i,j)
         o(j,i)=a(j,i)
         end do
      end do  
end subroutine
