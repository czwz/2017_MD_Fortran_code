program fcc_mov
use g
implicit none
character(len=9)::      FILE1,FILE2
integer::               T
real*8,allocatable::    xc(:),yc(:),zc(:),vcx(:),vcy(:),vcz(:)

        CALL IMPORTDATA
        CALL EQFCCFILE
        allocate(xc(atom),yc(atom),zc(atom),vcx(atom),vcy(atom),vcz(atom))

do T=40,160,40

        CALL INITIALIZE(T,xc,yc,zc,vcx,vcy,vcz)
        write(FILE1,'("Traj_",i3,"K")') T
        write(FILE2,'("Engy_",i3,"K")') T
        open(unit=1,file=FILE1,status='unknown')
        open(unit=2,file=FILE2,status='unknown')
        CALL VELOCITY_VERLET(xc,yc,zc,vcx,vcy,vcz,T)
        close(1,status='keep')
        close(2,status='keep')

enddo

stop
end program
