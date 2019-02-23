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

subroutine INITIALIZE(T,xc,yc,zc,vcx,vcy,vcz)
real*8,intent(inout):: xc(:),yc(:),zc(:),vcx(:),vcy(:),vcz(:)
real*8::               A,B,C,v1(atom),v2(atom),v3(atom),v_factor,Ek
integer::              i
integer,intent(in)::   T
mass=39.948                                              !amu
timestep=0.001                                           !ps 
gg=3*atom                                                !For Hoover Thermostat in real sampling      
open(unit=1,file='fcc',status='unknown')
        CALL random_seed()
do i=1,atom
        read(1,*) xc(i),yc(i),zc(i)
        CALL random_number(A)
        CALL random_number(B)
        CALL random_number(C)
        vcx(i)=2*A-1
        vcy(i)=2*B-1
        vcz(i)=2*C-1
        Ek=Ek+0.5*mass*(vcx(i)**2+vcy(i)**2+vcz(i)**2)
enddo
        v_factor=(0.5*gg*kb*real(T)*6.02*1000/Ek)**0.5
do i=1,atom
        vcx(i)=vcx(i)*v_factor*(0.01)
        vcy(i)=vcy(i)*v_factor*(0.01)
        vcz(i)=vcz(i)*v_factor*(0.01)
enddo
end subroutine

subroutine VELOCITY_VERLET(xc,yc,zc,vcx,vcy,vcz,T)
real*8,intent(inout)::  xc(:),yc(:),zc(:),vcx(:),vcy(:),vcz(:)
real*8,allocatable::    fx(:),fy(:),fz(:),xf(:),yf(:),zf(:),vhalf_x(:),vhalf_y(:),vhalf_z(:)
real*8::                r,factor,Ek,Untot,T0,Ti
integer::               i,j,k,l,iter
integer,intent(in)::    T
iter=5000
factor=eV*6.02*1000                                                         !anstromg/ps2
cutoff=2.5*ae                                                               !anstromg
T0=real(T)
Ti=real(T)
        allocate(fx(atom),fy(atom),fz(atom),xf(atom),yf(atom),zf(atom),vhalf_x(atom),vhalf_y(atom),vhalf_z(atom))
do i=1,iter
        Untot=0
        Ek=0
        write(1,*) atom
        write(1,*) "  "
                CALL force_calculator(xc,yc,zc,fx,fy,fz)
                CALL HOOVER_STAT(vcx,vcy,vcz,fx,fy,fz,T0)
        do j=1,atom
                vhalf_x(j)=vcx(j)+0.5*factor*(fx(j)/mass)*(timestep)
                vhalf_y(j)=vcy(j)+0.5*factor*(fy(j)/mass)*(timestep)
                vhalf_z(j)=vcz(j)+0.5*factor*(fz(j)/mass)*(timestep)
                xf(j)=xc(j)+vhalf_x(j)*timestep
                yf(j)=yc(j)+vhalf_y(j)*timestep
                zf(j)=zc(j)+vhalf_z(j)*timestep
                CALL MIRROR(xf(j),yf(j),zf(j))
                Ek=Ek+0.5*mass*(vcx(j)**2+vcy(j)**2+vcz(j)**2)*((10**(-3.))/(6.02*eV))
                write(1,*) "Ar",xc(j),yc(j),zc(j)
                do l=(j+1),atom
                        r=(pbc((xc(j)-xc(l)),ae)**2+pbc((yc(j)-yc(l)),ae)**2+pbc((zc(j)-zc(l)),ae)**2)**0.5
                                if (r.LE.cutoff) then
                                         Untot=Untot+4*epsilon*((sigma/r)**12-(sigma/r)**6)
                                endif   
                enddo
       enddo   
                100 format(I5," ",f10.5," ",f10.5," ",f10.5," ",f10.5," ",f20.10)
                T0=(Ek)*(eV/(kb*(3*real(atom))))*10000*2
                write(2,100) i,Untot,Ek,Untot+Ek,T0  
                CALL force_calculator(xf,yf,zf,fx,fy,fz)
        do j=1,atom
                vcx(j)=vhalf_x(j)+0.5*factor*(fx(j)/mass)*timestep
                vcy(j)=vhalf_y(j)+0.5*factor*(fy(j)/mass)*timestep
                vcz(j)=vhalf_z(j)+0.5*factor*(fz(j)/mass)*timestep
        enddo
                CALL HOOVER_STAT(vcx,vcy,vcz,fx,fy,fz,T0)
        do j=1,atom
                vcx(j)=vhalf_x(j)+0.5*factor*(fx(j)/mass)*timestep
                vcy(j)=vhalf_y(j)+0.5*factor*(fy(j)/mass)*timestep
                vcz(j)=vhalf_z(j)+0.5*factor*(fz(j)/mass)*timestep
                xc(j)=xf(j)
                yc(j)=yf(j)
                zc(j)=zf(j)
        enddo
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

subroutine HOOVER_STAT(vcx,vcy,vcz,fx,fy,fz,T)
real*8,intent(in)::    vcx(:),vcy(:),vcz(:)
real*8,intent(inout):: fx(:),fy(:),fz(:),T
real*8::               px(atom),py(atom),pz(atom),p,f,sum1,sum2,alpha
!integer,intent(in)::   T
integer::              i
sum1=0
sum2=0
do i=1,atom
       px(i)=mass*vcx(i)
       py(i)=mass*vcy(i)
       pz(i)=mass*vcz(i)
       sum1=sum1+(fx(i)*px(i)+fy(i)*py(i)+fz(i)*pz(i))/(2*mass)
       sum2=sum2+(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))/(2*mass)
enddo
       alpha=(sum1/sum2)
       !alpha=10000.*(sum1*eV/(0.5*(3*real(atom)-1)*kb*T))      !1/ps
do i=1,atom
       !fx(i)=fx(i)-alpha*px(i)/(6.02*eV*1000.)               !eV/anstromg
       !fy(i)=fy(i)-alpha*py(i)/(6.02*eV*1000.)
       !fz(i)=fz(i)-alpha*pz(i)/(6.02*eV*1000.)
       fx(i)=fx(i)-alpha*px(i)
       fy(i)=fy(i)-alpha*py(i)
       fz(i)=fz(i)-alpha*pz(i)
enddo
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
