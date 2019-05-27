program dipole_moment
  implicit none
  integer::nframe,nskip,x,nmols_per_stack,nstacks,iframe
  integer::nmols,natoms,tot_atoms
  real*8,allocatable::r(:,:),charge(:),dipole(:)
  character*25::file1,dumpfile,filename
  integer::i,j,k,n,l,dummy1,dummy2,frequency
  real*8::timestep,dt,polarization,DP,DP_xy,DP_z
  logical::exist
  real(kind=8)::box1,box2,rij(3),box_len,diff(3),avec,bvec,cvec,alpha,beta,gamma
  real(kind=8)::srij(3),h_mat(3,3),h_inv(3,3),vol
  real(kind=8)::xl,xh,yl,yh,zl,zh,xl_bound,xh_bound,yl_bound,yh_bound,zl_bound,zh_bound 
  real(kind=8)::xy,xz,yz,lx,ly,lz


  interface
  subroutine h_matrix(avec,bvec,cvec,alpha,beta,gamma,h_mat,h_inv,vol)
   real(kind=8)::avec,bvec,cvec,alpha,beta,gamma
   real(kind=8)::h_mat(3,3),h_inv(3,3),vol
  end subroutine h_matrix
 end interface

 
   call getarg(1,file1)
   open(2,file=file1,action="read")
   read(2,*)nmols
   read(2,*)natoms
   read(2,*)nmols_per_stack,nstacks
   read(2,*)timestep,frequency
   read(2,*)nframe,nskip
   read(2,*)dumpfile
   read(2,*)x
 
   ! convertion from fs to ns
   dt=real(frequency)*timestep*0.001d0*0.001d0
   tot_atoms=(nmols*natoms)
   
   do i=1,nstacks
       write(filename,'("DP",I1,".xvg")')i
       inquire(file=filename, exist=exist)
       if (exist) then
            open(40+i,file=filename,status="old",position="append",action="write")      ! bulk dipolemoment in a stack
       else
            open(40+i,file=filename,status="new",action="write")      ! bulk dipolemoment in a stack
            write(40+i,*)"@ title ""Total Dipole Moment "" "
            write(40+i,*)"@ xaxis label ""Time (ns)"""
            write(40+i,*)"@ yaxis label ""Dipole moment (D)"""
            write(40+i,*)"@ TYPE xy"
            write(40+i,*)"@ view 0.15, 0.15, 0.75, 0.85"
            write(40+i,*)"@ legend on"
            write(40+i,*)"@ legend box on"
            write(40+i,*)"@ legend loctype view"
            write(40+i,*)"@ legend 0.78, 0.8"
            write(40+i,*)"@ legend length 2"
            write(40+i,"(A4,I1,A16,I1,A3)")"@ s",i-1," legend ""Stack\s",i,"\N"""
       end if
       write(filename,'("DP_z",I1,".xvg")')i
       inquire(file=filename, exist=exist)
       if (exist) then
             open(100+i,file=filename,status="old",position="append",action='write')
       else
             open(100+i,file=filename,status="new",action='write')
             write(100+i,*)"@ title ""Total Dipole Moment along Stacking &
&Direction"" "
             write(100+i,*)"@ xaxis label ""Time (ns)"" "
             write(100+i,*)"@ yaxis label ""Dipole moment along Stacking Direction (D)"""
             write(100+i,*)"@ TYPE xy"
             write(100+i,*)"@ view 0.15, 0.15, 0.75, 0.85"
             write(100+i,*)"@ legend on"
             write(100+i,*)"@ legend box on"
             write(100+i,*)"@ legend loctype view"
             write(100+i,*)"@ legend 0.78, 0.8"
             write(100+i,*)"@ legend length 2"
             write(100+i,"(A4,I1,A16,I1,A3)")"@ s",i-1," legend ""Stack\s",i,"\N"""
       endif

       write(filename,'("DP_xy",I1,".xvg")')i
       inquire(file=filename, exist=exist)
       if (exist) then
             open(200+i,file=filename,status="old",position="append",action='write')
       else
             open(200+i,file=filename,status="new",action='write')
             write(200+i,*)"@ title ""Total Dipole Moment along Stacking &
&Direction"" "
             write(200+i,*)"@ xaxis label ""Time (ns)"" "
             write(200+i,*)"@ yaxis label ""Dipole moment in XY (D)"""
             write(200+i,*)"@ TYPE xy"
             write(200+i,*)"@ view 0.15, 0.15, 0.75, 0.85"
             write(200+i,*)"@ legend on"
             write(200+i,*)"@ legend box on"
             write(200+i,*)"@ legend loctype view"
             write(200+i,*)"@ legend 0.78, 0.8"
             write(200+i,*)"@ legend length 2"
             write(200+i,"(A4,I1,A16,I1,A3)")"@ s",i-1," legend ""Stack\s",i,"\N"""
       endif
   enddo
  
   
   inquire(file="polarization.xvg", exist=exist)
   if (exist) then
       open(2000,file="polarization.xvg",status="old",position="append",action="write")      ! polarization
   else
       open(2000,file="polarization.xvg",status="new",action="write")      ! polarization 
       write(2000,*)"@ title ""Polarization"" "
       write(2000,*)"@ xaxis label ""Time (ns)"""
       write(2000,*)"@ yaxis label ""Polarization (\xm\f{}C/cm\S2\N)"""
       write(2000,*)"@ TYPE xy"
       write(2000,*)"@ view 0.15, 0.15, 0.75, 0.85"
   end if

   inquire(file="DP.xvg", exist=exist)
   if (exist) then
       open(3000,file="DP.xvg",status="old",position="append",action="write")      ! polarization
   else
       open(3000,file="DP.xvg",status="new",action="write")      ! polarization 
       write(3000,*)"@ title ""Dipole moment"" "
       write(3000,*)"@ xaxis label ""Time (ns)"""
       write(3000,*)"@ yaxis label ""Dipole moment (D)"""
       write(3000,*)"@ TYPE xy"
       write(3000,*)"@ view 0.15, 0.15, 0.75, 0.85"
       write(3000,*)"@ legend on"
       write(3000,*)"@ legend box on"
       write(3000,*)"@ legend loctype view"
       write(3000,*)"@ legend 0.78, 0.8"
       write(3000,*)"@ legend length 2"
       write(3000,*)"@ s0 legend ""DP"""
       write(3000,*)"@ s1 legend ""DP\sxy\N"""
       write(3000,*)"@ s2 legend ""DP\sz\N"""
   end if


   allocate(charge(natoms))
   allocate(r(tot_atoms,3),dipole(3))

   do i=1,natoms
      read(2,*)charge(i)
   enddo
   
   open(8,file=dumpfile,action='read')


   do i=1,nskip+x
       do j=1,tot_atoms+9
           read(8,*)
       enddo
   enddo

   nframe=nframe-nskip-x

   
      ! calculating the dipole in the stack
   
  do iframe=1,nframe
   polarization=0.0d0;DP=0.d0;DP_xy=0.0d0;DP_z=0.0d0
   write(*,"(A15,3x,I5,A1,I5)")"Frame Progress:",iframe,"/",nframe
   ! Reading the dumpfile  ........1......................................................................
       do k=1,5
              read(8,*)
       enddo
       read(8,*)xl_bound,xh_bound,xy
       read(8,*)yl_bound,yh_bound,xz
       read(8,*)zl_bound,zh_bound,yz
       read(8,*)


       xl=xl_bound-min(0.0,xy,xz,xy+xz)
       xh=xh_bound-max(0.0,xy,xz,xy+xz) 
       yl=yl_bound-min(0.0,yz)
       yh=yh_bound-max(0.0,yz)
       zl=zl_bound
       zh=zh_bound
       lx=xh-xl
       ly=yh-yl
       lz=zh-zl
       avec=lx
       bvec=sqrt(ly**2+xy**2)
       cvec=sqrt(lz**2+xz**2+yz**2)
       alpha=acos((xy*xz+ly*yz)/(bvec*cvec))
       beta=acos(xz/cvec)
       gamma=acos(xy/bvec)
       call  h_matrix(avec,bvec,cvec,alpha,beta,gamma,h_mat,h_inv,vol)

       do j=1,tot_atoms
           read(8,*)n,dummy1,dummy2,r(n,:)
       enddo
! Calculating the dipole moment
        do i=1,nstacks
           dipole=0.0d0;l=0
           do j=1,nmols_per_stack
               do k=1,natoms
                  l=l+1
                  n=((i-1)*natoms*nmols_per_stack)+l
                 ! write(*,*)"nstack,nmols_per_stack,natoms",i,j,k,n
                  dipole(:)=dipole(:)+r(n,:)*charge(k)
               enddo
           enddo  ! End of stacks loop
           dipole(:)=dipole(:)/0.20819434d0
           write(100+i,1000)((iframe+x+nskip-1)*dt),dipole(3)
           write(200+i,1000)((iframe+x+nskip-1)*dt),dsqrt(dipole(2)**2+dipole(1)**2)
           write(40+i,1000)((iframe+x+nskip-1)*dt),dsqrt(dot_product(dipole(:),dipole(:)))
           polarization=polarization+dsqrt(dot_product(dipole(:),dipole(:)))
           DP=dsqrt(dot_product(dipole(:),dipole(:)))+DP
           DP_xy=dsqrt(dipole(2)**2+dipole(1)**2)+DP_xy
           DP_z=dipole(3)+DP_z
        enddo ! End of all stacks
        polarization=polarization/vol
        write(2000,1000)((iframe+x+nskip-1)*dt),polarization*333.563d0
        write(3000,10000)((iframe+x+nskip-1)*dt),DP/0.20819434d0,DP_xy/0.20819434d0,DP_z/0.20819434d0
   enddo   
1000 format(F15.7,3x,F15.7)
10000 format(F15.7,3x,3(F15.7,3x))
end program dipole_moment


subroutine h_matrix(avec,bvec,cvec,alpha,beta,gamma,h_mat,h_inv,vol)
 implicit none
 real(kind=8)::avec,bvec,cvec,alpha,beta,gamma,h_mat(3,3),h_inv(3,3),vol
 h_mat(1,1)=avec
 h_mat(1,2)=0.d0
 h_mat(1,3)=0.d0
 h_mat(2,1)=bvec*cos(gamma)
 h_mat(2,2)=bvec*sin(gamma)
 h_mat(2,3)=0.d0
 h_mat(3,1)=cvec*cos(beta)
 h_mat(3,2)=cvec*(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)
 h_mat(3,3)=cvec*sqrt(((sin(beta))**2)-(h_mat(3,2)/cvec)**2)
 h_mat=transpose(h_mat)
 vol=h_mat(1,1)*(h_mat(2,2)*h_mat(3,3)-h_mat(2,3)*h_mat(3,2)) &
-h_mat(1,2)*(h_mat(2,1)*h_mat(3,3)-h_mat(2,3)*h_mat(3,1)) &
       +h_mat(1,3)*(h_mat(2,1)*h_mat(3,2)-h_mat(3,1)*h_mat(2,2)) 
 h_inv(1,1)=(h_mat(2,2)*h_mat(3,3)-h_mat(2,3)*h_mat(3,2))/(vol)
 h_inv(1,2)=(h_mat(1,3)*h_mat(3,2)-h_mat(1,2)*h_mat(3,3))/(vol)
 h_inv(1,3)=(h_mat(1,2)*h_mat(2,3)-h_mat(1,3)*h_mat(2,2))/(vol)
 h_inv(2,1)=(h_mat(2,3)*h_mat(3,1)-h_mat(2,1)*h_mat(3,3))/(vol)
 h_inv(2,2)=(h_mat(1,1)*h_mat(3,3)-h_mat(1,3)*h_mat(3,1))/(vol)
 h_inv(2,3)=(h_mat(1,3)*h_mat(2,1)-h_mat(1,1)*h_mat(2,3))/(vol)
 h_inv(3,1)=(h_mat(2,1)*h_mat(3,2)-h_mat(2,2)*h_mat(3,1))/(vol)
 h_inv(3,2)=(h_mat(1,2)*h_mat(3,1)-h_mat(1,1)*h_mat(3,2))/(vol)
 h_inv(3,3)=(h_mat(1,1)*h_mat(2,2)-h_mat(2,1)*h_mat(1,2))/(vol)
 return
end subroutine h_matrix
