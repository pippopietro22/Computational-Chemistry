program Second_Quantization
    implicit none

!   pos=position to occupy - i=counter- n=number of eletrons counted - segno=factor filling order
    integer :: pos,i,n,segno

!   k=binary vector spin-orbital - s=confrontation vector - mask=electron counting vector
    integer :: k,s,mask

!   ciclo=parameter for iteration filling
    character :: ciclo
!
    ciclo='y'     !until ciclo='y' the iteration continues
    segno=1       

!   Insert starting vector (decimal base)
    write(*,*)'Insert starting vector in decimal base:'
    read(*,*)k
    write(*,'(b64)')k
    write(*,*)''
!    
    do while (ciclo .eq. 'y')
        write(*,*)'Spin-Orbitale to fill?'
        write(*,*)'(Insert a number grater than 100 to exit the program)'
        read(*,*)pos
        write(*,*)''
!
        if (pos .le. 25)then             !You can fill until 25th position

!           Checks if the spin-orbital is already filled
            s=ibset(0,pos-1)
            if(and(s,k) .ne. s)then

!               Checks how many spin-orbitals are already filled
                mask=0
                do i=1,pos-2
                    mask=ibset(mask,i)
                end do
            
!               Calculating fase of the function
                n=popcnt(and(mask,k))
                segno=segno*((-1)**n)

!               Finally, inserts the electron in the wanted position
                k=ibset(k,pos-1)
                write(*,'(b64)')k
                write(*,*)'factor',segno
                write(*,*)'in decimal base is:',k
                write(*,*)''
            else

!               if the position is already filled, the program doesn't stop, it just tells you and keep running
                write(*,*)'Position already occupied.'
                write(*,*)''
            end if
        else

!           If you choose a position superior to 25, the program ends
            ciclo='n'
        end if
    end do

!   prints the results both in binary and decimal base with his fase
    open(11,file='Function spin-orbital.txt')
        write(11,'(b64)')k
        write(11,*)'factor',segno
        write(11,*)'in decimal base is',k
    close(11)
    
end program
