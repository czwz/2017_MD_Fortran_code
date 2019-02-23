module g
implicit none
integer,public    ::line,atom,gg
real*8,public     ::sigma,epsilon,cutoff,ae,mass,timestep,kb,eV
contains

subroutine IMPORTDATA
sigma=3.40
epsilon=0.01
kb=1.3806485279
eV=1.60217662
line=5
end subroutine 

subroutine fcc_generator(a)
real*8,intent(in):: a
integer:: i,j,k
        open(unit=2,file='fcc',status='unknown')
        do i=0,(line-1)
                do j=0,(line-1)
                        do k=0,(line-1)
                                write(2,*) (0+i)*a,(0+j)*a,(0+k)*a
                                write(2,*) (0.5+i)*a,(0.5+j)*a,(0+k)*a
                                write(2,*) (0+i)*a,(0.5+j)*a,(0.5+k)*a
                                write(2,*) (0.5+i)*a,(0+j)*a,(0.5+k)*a
                        enddo
                enddo
        enddo
        close(2,status='keep')
end subroutine

subroutine EQFCCFILE
real*8,allocatable::   xx(:),yy(:),zz(:)
real*8,allocatable::   rn(:),Un(:)
real*8::               r,a,space,U
integer::              i,j,k,step
space=(1/1000.)
atom=(line**3)*4
        allocate(xx(atom),yy(atom),zz(atom))
        open(unit=3,file='Un-r',status='unknown')
        step=5000
        do k=1,step
                a=3.+space*k
                cutoff=2.5*a

                CALL fcc_generator(a)
                open(unit=2,file='fcc',status='unknown')
                do i=1,atom 
                        read(2,*) xx(i),yy(i),zz(i)
                enddo
                close(2,status='keep')
                U=0
                do i=1,atom
                        do j=(i+1),atom
                                r=(pbc((xx(i)-xx(j)),a)**2+pbc((yy(i)-yy(j)),a)**2+pbc((zz(i)-zz(j)),a)**2)**0.5
                                if (r.LE.cutoff) then
                                        U=U+4*epsilon*((sigma/r)**12-(sigma/r)**6)
                                endif
                        enddo
                enddo
            write(3,*) a,U
        enddo
        close(3,status='keep')

        open(unit=3,file='Un-r',status='unknown')
        allocate(rn(step),Un(step))
        do i=1,step
                read(3,*) rn(i),Un(i)
        enddo
        U=minval(Un)
        do i=1,step
                if (Un(i).EQ.U) then
                        ae=rn(i)
                        CALL fcc_generator(ae)
                endif
        enddo
        close(3,status='keep')
end subroutine

subroutine INITIALIZE(T,xc,yc,zc,xp,yp,zp)
real*8,intent(inout):: xc(:),yc(:),zc(:),xp(:),yp(:),zp(:)
real*8::               A,B,C,v1(atom),v2(atom),v3(atom),v_factor,Ek
integer::              i
integer,intent(in)::   T
mass=39.948                                              !amu
timestep=0.001                                           !ps 
gg=3*atom-1                                              !For Hoover Thermostat in real sampling      
        OPEN(unit=1,file='fcc',status='unknown')
do i=1,atom
        read(1,*) xc(i),yc(i),zc(i)                      !update current position rc
enddo
        CALL random_seed()
do i=1,atom
        CALL random_number(A)
        CALL random_number(B)
        CALL random_number(C)
        v1(i)=2*A-1
        v2(i)=2*B-1
        v3(i)=2*C-1
        Ek=Ek+0.5*mass*(v1(i)**2+v2(i)**2+v3(i)**2)
enddo
do i=1,atom
        v_factor=(0.5*gg*kb*real(T)*6.02*1000/Ek)**0.5
        xp(i)=-v1(i)*v_factor*(0.01)*timestep+xc(i)
        yp(i)=-v2(i)*v_factor*(0.01)*timestep+yc(i)
        zp(i)=-v3(i)*v_factor*(0.01)*timestep+zc(i)
enddo
end subroutine

subroutine VERLET(xc,yc,zc,xp,yp,zp,T)
real*8,intent(inout)::  xc(:),yc(:),zc(:),xp(:),yp(:),zp(:)
real*8,allocatable::    fx(:),fy(:),fz(:),xf(:),yf(:),zf(:),Untot(:),Ek(:),v(:)
real*8::                r,factor,vx,vy,vz,v_factor,E,T0,Ti,alpha
integer::               i,j,k,l,iter
integer,intent(inout):: T
iter=5000
factor=eV*6.02*1000                                                         !anstromg/ps2
cutoff=2.5*ae                                                               !anstromg
T0=real(T)
Ti=real(T)
        allocate(fx(atom),fy(atom),fz(atom),xf(atom),yf(atom),zf(atom),Untot(iter),Ek(iter),v(atom))
do i=1,iter
        Untot(i)=0
        Ek(i)=0
        write(1,*) atom
        write(1,*) "  "
                CALL force_calculator(xc,yc,zc,fx,fy,fz)
        do j=1,atom
                xf(j)=2*xc(j)-xp(j)+factor*(fx(j)/mass)*(timestep)**(2)     !update future position rf
                yf(j)=2*yc(j)-yp(j)+factor*(fy(j)/mass)*(timestep)**(2)
                zf(j)=2*zc(j)-zp(j)+factor*(fz(j)/mass)*(timestep)**(2)
                CALL MIRROR(xf(j),yf(j),zf(j))
        enddo
        do j=1,atom
                vx=(pbc(xf(j)-xc(j),ae))/(timestep)                       !anstromg/ps
                vy=(pbc(yf(j)-yc(j),ae))/(timestep)
                vz=(pbc(zf(j)-zc(j),ae))/(timestep)
                Ek(i)=Ek(i)+0.5*mass*(vx**2+vy**2+vz**2)*(10**(-3.)/(6.02*eV))              !eV
        enddo
                alpha=((0.5*(3*real(atom))*kb*real(Ti)/(10000*Ek(i)*eV))**0.5)
        do j=1,atom
                write(1,*) "Ar",xc(j),yc(j),zc(j)
                xp(j)=xf(j)-(pbc(xf(j)-xc(j),ae)*alpha)
                yp(j)=yf(j)-(pbc(yf(j)-yc(j),ae)*alpha)
                zp(j)=zf(j)-(pbc(zf(j)-zc(j),ae)*alpha)
                xc(j)=xf(j)
                yc(j)=yf(j)
                zc(j)=zf(j)
        enddo
                do k=1,atom
                        do l=(k+1),atom
                                r=(pbc((xp(k)-xp(l)),ae)**2+pbc((yp(k)-yp(l)),ae)**2+pbc((zp(k)-zp(l)),ae)**2)**0.5
                                if (r.LE.cutoff) then
                                        Untot(i)=Untot(i)+4*epsilon*((sigma/r)**12-(sigma/r)**6)
                                endif
                        enddo
                enddo
                100 format(I5," ",f10.5," ",f10.5," ",f10.5," ",f10.5," ",f20.10)
                T0=(Ek(i))*(eV/(kb*(3*real(atom)-1)))*10000*2
                write(2,100) i,Untot(i),Ek(i),Untot(i)+Ek(i),T0                                              !eV, eV, eV, K
enddo
end subroutine

subroutine MIRROR(xx,yy,zz)
real*8,intent(inout):: xx,yy,zz
        if (xx.GT.ae*line) then
                xx=xx-ae*line
        elseif (xx.LT.0) then
                xx=xx+ae*line
        endif

        if (yy.GT.ae*line) then
                yy=yy-ae*line
        elseif (yy.LT.0) then
                yy=yy+ae*line
        endif

        if (zz.GT.ae*line) then
                zz=zz-ae*line
        elseif (zz.LT.0) then
                zz=zz+ae*line
        endif
end subroutine

subroutine force_calculator(x,y,z,fx,fy,fz)
real*8,intent(inout):: x(:),y(:),z(:),fx(:),fy(:),fz(:)
real*8:: r,rx,ry,rz,a1,a2,a3
integer:: i,j
do i=1,atom
        fx(i)=0
        fy(i)=0
        fz(i)=0
        do j=1,atom
                if (i.NE.j) then
                        !rx=(x(j)-x(i))                 !Cluster
                        !ry=(y(j)-y(i))
                        !rz=(z(j)-z(i))
                        rx=pbc(x(j)-x(i),ae)
                        ry=pbc(y(j)-y(i),ae)
                        rz=pbc(z(j)-z(i),ae)
                        r=(rx**2+ry**2+rz**2)**0.5
                                if (r.LE.cutoff) then
                                        fx(i)=fx(i)-(48*epsilon/sigma)*((sigma/r)**13-0.5*(sigma/r)**7)*(rx/r)
                                        fy(i)=fy(i)-(48*epsilon/sigma)*((sigma/r)**13-0.5*(sigma/r)**7)*(ry/r)
                                        fz(i)=fz(i)-(48*epsilon/sigma)*((sigma/r)**13-0.5*(sigma/r)**7)*(rz/r)
                                endif
                endif
        enddo
enddo
end subroutine

function pbc(r,a)
real*8,intent(in):: r,a
real*8:: pbc,l
l=((line)*a)
        if ((r.GT.0).AND.(r**2.GT.(0.5*l)**2)) then
                pbc=(r-l)
        elseif ((r.LT.0).AND.(r**2.GT.(0.5*l)**2)) then 
                pbc=(r+l)
        elseif ((r**2).LT.(0.5*l)**2) then
                pbc=(r)
        endif
return
end function

end module
