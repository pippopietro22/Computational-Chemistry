Program Cholesky
        implicit none
        Real*8 :: somma
        Integer :: n,i,j,k
        Real*8, dimension(:,:), allocatable :: a,l,g
!
!       Matrix Dimension
!
        write(*,*)'Insert (square) matrix dimension: '
        read(*,*)n
!
        allocate (a(n,n),l(n,n),g(n,n))
!
!       Matrix Generation (You can insert in this section your own matrix to decompose)
!
        do i=1,n
                a(i,i)=real(i,8) + 1.d0
                do j=1,i-1
                        
                        a(j,i)=1.d0/(real(j,8)+real(i,8))
                        a(i,j)=a(j,i)
                end do
        end do    
!
!       Cholensky Algorithm
!
        do i=1,n 
                do j=1,i
                        if (i .eq. j) then
                                somma=0.d0
                                do k=1,j-1
                                        somma = somma + l(i,k)**2
                                end do 
                                l(j,j)=sqrt(a(j,j)-somma)
                        else
                                somma=0.d0
                                do k=1,j-1
                                        somma=somma+l(i,k)*l(j,k)
                                end do
                                l(i,j)=(a(i,j) - somma)/l(j,j)
                        end if
                end do 
        end do
!
!       Transpose Matrix
!
        do i=1,n
                do j=1,n 
                        g(j,i)=l(i,j)
                end do 
        end do 
!
!       Writing resulting Matrix
!
        do i=1,n
                write(*,*)(a(i,j), j=1,n)
        end do
        write(*,*)''
        do i=1,n
                write(*,*)(l(i,j), j=1,n)
        end do
end program Cholesky        
