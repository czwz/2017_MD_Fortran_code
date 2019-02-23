program fcc_mov
use g
implicit none
character(len=9)::      FILE1,FILE2
integer::               T
real*8,allocatable::    xc(:),yc(:),zc(:),xp(:),yp(:),zp(:)

        CALL IMPORTDATA
        CALL EQFCCFILE
        allocate(xc(atom),yc(atom),zc(atom),xp(atom),yp(atom),zp(atom))

T=80

        CALL INITIALIZE(T,xc,yc,zc,xp,yp,zp)
        write(FILE1,'("Traj_",i3,"K")') T
        write(FILE2,'("Engy_",i3,"K")') T
        open(unit=1,file=FILE1,status='unknown')
        open(unit=2,file=FILE2,status='unknown')
        CALL VERLET(xc,yc,zc,xp,yp,zp,T)
        close(1,status='keep')
        close(2,status='keep')



stop
end program
