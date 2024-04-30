program CG
        implicit none
        integer :: i,j,n,g
        Real*8 :: alfa, beta, conv, tau, gamma, gammaold
        Real*8, dimension(:,:), allocatable :: a
        Real*8, dimension(:), allocatable :: x,x1,b,r,r1,h
        !
        write(*,*)'Insert System dimension: '
        read(*,*)n
        allocate (a(n,n),x(n),x1(n),b(n),r(n),r1(n),h(n))
        !
        !preamble
        !
        call Pream (i,j,n,conv,a,x,b,r,g,tau)
        !
        !Iniztializing Conjugate Gradient
        !
        call dgemv('N',n,n,1.d0,a,n,x,1,0.d0,r1,1)
        !
        do i=1,n
                r1(i)=r1(i)-b(i)
        end do
        !
        beta=0.d0
        !
        gammaold=dot_product(r1,r1)
        !
        do while (tau .gt. conv)
                g=g+1
                !
                !Alpha term creation
                !
                call dgemv('N',n,n,1.d0,a,n,r1,1,0.d0,h,1)
                !
                gamma=dot_product(r1,h)
                !
                alfa=-(gammaold/gamma)
                !
                !x1 vector creation
                !
                do i=1,n
                        x1(i)=x(i) + alfa*(r1(i)+beta*r(i))
                end do
                !
                !saving old value for r1 in r
                !
                do i=1,n
                        r(i)=r1(i)
                end do
                !
                !Calculating r1
                !
                call dgemv('N',n,n,1.d0,a,n,x1,1,0.d0,h,1)
                !
                do i=1,n
                        r1(i)=h(i)-b(i)
                end do
                !
                !Calculating residual
                !
                do i=1,n
                        tau=r1(i)**2
                end do
                tau=sqrt(tau)
                !
                !Beta term creation
                !
                gamma=dot_product(r1,r1)
                beta=gamma/gammaold
                !
                !Saving variables
                !
                gammaold=gamma
                !
                do i=1,n
                        x(i)=x1(i)
                end do
                !
                if (g .eq. n) then
                        tau = 0.0d0
                end if
        end do
        write(*,*)'Number of iterations: ', g
end program
!
!preambolo
!
subroutine Pream (i,j,n,conv,a,x,b,r,g,tau)
        implicit none 
        integer :: i,j,n,g
        Real*8 :: conv,tau
        Real*8 :: a (n,n),x(n),b(n),r(n)
        !
        !fattore convergenza
        !
        write(*,*)'Insert convergence factor: '
        read(*,*)conv
        !
        !creazione matrice a
        !
        call SuperMatrixMaker (i,j,n,a)
        !
        !creazione termine noto b
        !
        do i=1,n
                b(i)=1/real(i,8) + 3.d0 
        end do
        !
        !creazione vettore x
        !
        do i=1,n
                x(i)=0.d0
        end do
        !
        !Inizializzazione variabili
        !
        do i=1,n
                r(i)=0.d0
        end do
        g=0
        tau=10*conv
end subroutine
!
!creazione matrice a
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
