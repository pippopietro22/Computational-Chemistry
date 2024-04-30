Program Davidson
    Implicit none
!   Indices/Cycles counters - g=counts number of iterations
!   controllo=allows to exit cycle when Nomra_residui is lower than tau
    Integer :: i,j,k,g,controllo

!   n=system dimension - giri=tells wich vector stock we are saving in matrix V - dimH=reduced matrix (H) dimension
    Integer :: n,giri,nvec,Vcolonne,dimH           

!   iter_max=maximum number of iterations    
    Integer, parameter :: Iter_max=100    
    
!   lwork and info are variables for lapack subroutines  
    Integer :: lwork,info     
    
!   tau=convergence - Residuo=parameter for residuals  
    real*8 :: tau,Residuo,lw(1), xx(1)   
    
!   eigenvalue=for reduced matrix - diagonal=matrix A diagonal - Norma_residui =  norm for vector "Residui"
    real*8, dimension(:), allocatable :: eigvalue,diagonal, Norma_residui

!   other vectors for lapak subroutines 
    real*8, dimension(:), allocatable :: taux,work  

!   A=Initial Matrix - V=Canonic Base - W=AV    
    Real*8, dimension(:,:), allocatable :: A,V,traspostaV,W    

!   r=residuals - Ritzvec - H= A projection on V - y=new vectors  
    Real*8, dimension(:,:), allocatable :: r,Ritzvec,H,U,y  

!   proiezY=Y projection on V - ComponentiYV=Y components during orthogonalization to V    
    Real*8, dimension(:,:), allocatable :: proiezY,ComponentiYV 
!
    write(*,*)'Insert dimension:'
    read(*,*) n
    if(n .lt. 100)then
        write(*,*)'This program is not capable do work with such a small system:'
        stop 10
    end if
!
    write(*,*)'Insert number of eigenvectors and eigen value to find:'
    read(*,*) nvec
    if(nvec .gt. int(n/10))then
        write(*,*)'TOO MANY requested.'
        stop 11
    end if

!   Vcolonna=Effective number of V columns (orthonormal base)
    Vcolonne=20*nvec        
    tau=1E-8
!
    lw=0.d0
    call dsyev('v','l',Vcolonne,xx,Vcolonne,xx,lw,-1,info)
    lwork = int(lw(1))
!
    allocate(A(n,n),V(n,Vcolonne),traspostaV(Vcolonne,n),W(n,Vcolonne),H(Vcolonne,Vcolonne),U(Vcolonne,Vcolonne),&
            Ritzvec(n,nvec),r(n,nvec),y(n,nvec)) 
    allocate(diagonal(n),eigvalue(Vcolonne),Norma_residui(nvec))       

!   Lapack Parameters        
    allocate(taux(n),work(lwork),proiezY(Vcolonne,nvec),ComponentiYV(n,nvec))  
!
!   Creating diagonal dominating matrix
!
    call MatrixMaker (n,A,diagonal)                           
!
!   Initializing subspace V
!
    v=0.d0
    do i=1,nvec
        v(i,i)=1.d0
    end do
!
!   Initializing Davidson
! 
    giri=1             
    g=0                 
    do k=1,iter_max
        g=g+1
        dimH=giri*nvec
!
!       Creating matrix H
!
        call dgemm('n','n',n,dimH,n,1.d0,A,n,V,n,0.d0,W,n)
!
        call dgemm('t','n',dimH,dimH,n,1.d0,V,n,W,n,0.d0,H,Vcolonne)
!        
!       Diagonalizing matrix H
!
        U=H
        call dsyev('V','L',dimH,U,Vcolonne,eigvalue,work,lwork,info) 
        if (info.ne.0) then
            write(6,*) ' warning: dsyev failed, error ', info
            stop 12
        end if
!
        call dgemm('n','n',n,nvec,dimH,1.d0,V,n,U,Vcolonne,0.d0,Ritzvec,n)
!
        call dgemm('n','n',n,nvec,dimH,1.d0,W,n,U,Vcolonne,0.d0,r,n)
        do i=1,nvec
            r(:,i)= r(:,i) - eigvalue(i)*Ritzvec(:,i) 
        end do
!
        controllo=0
        do i=1,nvec
            Residuo=norm2(r(:,i))
            Norma_residui(i)=Residuo
!
            if(Residuo .gt. tau)then
                controllo=1
            end if
        end do
!
        if(controllo .eq. 0)then
            exit
        end if
!
!       residual preconditioning
!
        do j = 1, nvec
          do i = 1, n
            residuo = diagonal(i) - eigvalue(1)
            if (residuo .gt. 1.0d-5) then
              y(i,j) = r(i,j) / residuo
            else
              y(i,j) = r(i,j)
            end if
          end do
        end do
!
!       Y orthogonalization 

!         1) Removing from Y his projection on V base:
        call dgemm('t','n',dimH,nvec,n,1.0d0,v,n,y,n,0.0d0,proiezy,vcolonne)
!
!            Now lets do (y - v * proiezy):
!
        call dgemm('n','n',n,nvec,dimH,-1.0d0,v,n,proiezy,vcolonne,1.0d0,y,n)
!
!         2) Now we have to orthonormalize y:

!       QR Decomposition
        call dgeqrf(n,nvec,y,n,taux,work,lwork,info)
        if (info.ne.0) then
            write(6,*) ' warning: dgeqrf failed, error ', info
            stop 13
        end if
!
        call dorgqr(n,nvec,nvec,y,n,taux,work,lwork,info)
        if (info.ne.0) then
            write(6,*) ' Warrning: dorgqr failled, errore ', info
            stop 14
        end if
!
!         3) Because of malconditioning problems, we have to re-orthogonalize the y vectors
        call dgemm('t','n',dimH,nvec,n,1.0d0,v,n,y,n,0.0d0,proiezy,vcolonne)
        call dgemm('n','n',n,nvec,dimH,-1.0d0,v,n,proiezy,vcolonne,1.0d0,y,n)
!
!       Saving vectors
!
        if (giri .lt. 20)then
!
!           Saving new orthogonalized vectors
!
            do i=1,nvec
                j=dimH +i
                V(:,j)=y(:,i)
            end do
            giri=giri+1

        else            !if V matrix is completly filled, we have to empty it and start from the beginning with ritz vectors and the 
            giri=2      !new orthogonalized vectors
!
!           Saving ritz vectors in the V matrix
!
            v=0.d0
            do i=1,nvec
                V(:,i)=Ritzvec(:,i)
            end do
!
!           Saving new orthonormalized vectors
!
            do i=1,nvec
                j=nvec+i
                V(:,j)=y(:,i)/norm2(y(:,i))
            end do     
        end if  
    end do
end program   
!
!
!
!
!
!Creating matrix A
!
subroutine MatrixMaker (n,a,diagonal)
    implicit none
    Integer :: i, j, n
    Real*8 :: a(n,n),diagonal(n)
    do i=1,n
       a(i,i)=i+1.d0
       diagonal(i)=a(i,i)
       do j=1,i-1                      
          a(j,i)=1.d0/(real(j,8)+real(i,8))
          a(i,j)=a(j,i)
      end do
    end do  
 end subroutine