program Jacobi
        Implicit none
        Real*8, dimension(:,:), allocatable :: a, o
        Real*8, dimension(:), allocatable :: x, x1, b, d
        Real*8 :: somma2, dif, tau
        integer :: i,j,n,g
!
!       Insert Matrix Dimension
!
        write(*,*)'Insert Dimension: '
        read(*,*)n
!
!       Convergence Factor
!
        write(*,*)'Insert convergence factor: '
        read(*,*)tau
!
        allocate(a(n,n),o(n,n), x(n), x1(n), b(n), d(n))
!
!       Creating a positive defined matrix
!
        call SuperMatrixMaker (i,j,n,a)
!
!       Decomposition of matrix 'a' in d(diagonal) and o(off diagonal)
!
        call Matrix_dec (i,j,n,a,d,o)
!
!       Creating guess vector x
!
        do i=1,n
                x(i)=0.d0
        end do
!
!       Creating vector b (known terms of the linear system)
!
        do i=1,n
                b(i)=1/real(i,8) +3.d0
        end do
!
!       Jacobi iteration 
!
        call Jacobi_sub (tau, g,i,n,somma2,dif,x,x1,b,o,d)
!
        write(*,*)'Number of iterations: ', g
!
end program Jacobi
!
!subroutine: creates a positive defined matrix
!
subroutine SuperMatrixMaker (i,j,n,a)
        implicit none
        Integer :: i, j, n
        Real*8 :: a(n,n)
        do i=1,n
                a(i,i)=i+1.d0
                do j=1,i-1                      
                        a(j,i)=1.d0/(real(j,8)+real(i,8))
                        a(i,j)=a(j,i)
                end do
        end do  
end subroutine
!
!Subroutine: decomposition of matrix 'a' in d(diagonal) and o(off diagonal)
!
subroutine Matrix_dec (i,j,n,a,d,o)
        implicit none
        Integer :: i,j,n
        Real*8, dimension(n,n) :: a,o
        Real*8 :: d(n)
        do i=1,n
                d(i)=a(i,i)
        end do 
        !
        do i=1,n
                do j=1,n
                        if (i .eq. j)then
                                o(j,i)=0.d0
                        else
                                o(j,i)=a(j,i)
                        end if
                end do 
        end do 
end subroutine
!
!Jacobi
!
Subroutine Jacobi_sub (tau, g,i,n, somma2,dif,x,x1,b,o,d)
        implicit none
        Integer :: g,i,n
        Real*8 :: somma2, dif, tau
        Real*8, dimension(n) :: x,x1,b,y,d
        Real*8, dimension(n,n) :: o
        g=0
        somma2=10.d0
        do while (somma2 .gt. tau)
                g=g+1
                call dgemv('N',n,n,1.d0,o,n,x,1,0.d0,y,1)
                do i=1,n
                        x1(i)=(b(i)-y(i))/d(i)
                end do
!
!               calculation of the increment
!
                somma2=0.d0
                do i=1,n
                        dif=x(i)-x1(i)
                        somma2=somma2+dif**2
                end do
!
                somma2=sqrt(somma2)/n
!
                do i=1,n
                        x(i)=x1(i)
                end do    
        end do
end subroutine
