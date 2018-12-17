!C ****************************************************************************
! FILE: mpi_project.f
! DESCRIPTION:
!   2d laplace iteration on rectangular domain to demonstrate MPI 
! AUTHOR: Basma Tumi
!C LAST REVISED: 12/03/2018
!C ****************************************************************************
      program project
      include 'mpif.h'
      parameter (MASTER = 0)
      integer :: ims,ime,jms,jme,its,ite,jts,jte,ib,jb
! i=Row(y)index, j=col(x)index, iter=current iteration, 
!niter=#itereations,nstop=stopiteration.
      integer :: i,j,iter,niter, mxiter, nstop      
      parameter(mxiter=5000)
      integer :: numtasks, taskid, len, ierr 
      character(MPI_MAX_PROCESSOR_NAME) hostname
     double precision, dimension(:,:), pointer ::u,uold,f, res,u_exact
!integer :: M, N, ids, ide, jds, jde
      integer :: isize,jsize
! RIGHT=RightMsgtag, LEFT= LeftMsgtag, TOP=TopMsgTag, BOTTOM=BottomMsgTag.
       integer :: RIGHT, LEFT
      parameter (LEFT=100, RIGHT=101, TOP=110,BOTTOM=111)

      real  :: tol,h, d, err, gerr
      integer :: divideup
      tol=0.000001
 
!h=real((ite-its)/N)

      ids=1
      ide=1000
      jds=1
      jde=1100
      
        
      
      !integer numtasks, taskid, len, ierr
      !character(MPI_MAX_PROCESSOR_NAME) hostname
    
      
      call MPI_INIT(ierr)
      !who am I? What is my rank?
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      !find the largest integer N such that N^2 <= numtasks
      N=floor(sqrt(real(numtasks)))
      M = N
      print*, ' task ', taskid 
      print*, ' domain ',ids, ' : ',ide, ' jds ',':', jde
      print*, numtasks, ' proc grid ', M, ' by ', N
      ib=1+taskid/M
      ! M=5 taskid = 0: => jb = 1+(0)/5 = 1 + 0 = 1
      ! M=5 taskid = 4: => jb = 1+(4)/5 = 1 + 0 = 1 
      ! M=5 taskid = 5: => jb = 1+(5)/5 = 1 + 1 = 2
      jb=1+mod(taskid,M)

      ! have ids:ide jds:jde = domain dimensions
      ! compute isize, jsize = tile size in directions i and j
      ! number of points in i direction divided by number of 
      ! domains in i direction
      isize=divideup(ide-ids+1,M)
      jsize=divideup(jde-jds+1,N)
      !compute my ib, jb, its,ite,jts,jte
      its=ids+isize*(ib-1)
      jts=jds+jsize*(jb-1)
      !Example:  ib =1   its = 1
      ! ib =2   its = 1 + isize
      ite=min(ide,ids-1+isize*ib)
      jte=min(jde,jds-1+jsize*jb)
print*, ' task ', taskid
print*, ' block ',ib, jb 
print*, ' tile ', its, ':', ite, ':',jts, ':', jte


   !Set ims,ime,jms,jme
      ims=its-1
      ime=ite+1
      jms=jts-1
      jme=jte+1

h=real((ite-its)/N)
d=real((ide-ids)/N)

!initialize u,f
  do j=jms,jme
    do i=ims,ime
uold(i,j)=0
   end do
end do


!Exact Solution
  Do j=jms,jme
   do i=ims,ime
    f(i,j)=2*sin(i*h)*cos(j*h)
   end do
end do
!Saving the exact solutuion in one array u_exact

 Do j=jds,jde
   do i=ids,ide
    u_exact(i,j)=2*sin(i*d)*cos(j*d)
   end do
end do
print*, 'Exact Solution'
do i=ids,ide
    write(*,*) (u_exact(i,j), j=jds,jde)
end do



   
!Boundary Conditions


!left
if ( jb==1 )then 

do i=its,ite
  u(i,1)=sin(i*h)*cos(h)
end do
end if 
!Right
if (jb==N) then
do i=its,ite
   u(i,N)=sin(i*h)*cos(N*h)
end do
end if
!Top
if (ib==1) then
do j=jts,jte
    u(1,j)=sin(h)*cos(j*h)
end do
end if
!Bottom
if (ib==N) then
do j=jts,jte
   u(N,j)=sin(N*h)*cos(j*h)
end do
end if




!main Work

!all tasks get collective communication of niter from MAster: 
call MPI_Bcast(niter, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, IERR)
!WHERE count is 1 and datatype is mpi_integer.
!Do computations on block for niter


Do iter=1,niter !77
   Do j=jts,jte
      Do i=its,ite
          u(i,j)= 0.25*(h**2*f(i,j)+u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1))
                  
      End Do
   End Do
!Get localMaxAbsChange for this task and save old u.
 err=0
Do j=jts,jte
      Do i=its,ite
         
               err=max( abs(u(i,j)-uold(i,j)),err)
      uold(i,j)=u(i,j)
      End Do
   End Do
!MPI step: taskid "jb" sends right Boundry data to Right Neighbor "jb+1"
if (jb<N) then     !88
   call MPI_send(u(its, jte+1), &
& N, MPI_DOUBLE_PRECISION, Jb+1,RIGHT, MPI_COMM_WORLD, IERR )
!starting at address u(its,jte+1) for N Rows.
!MPI step: taskid "jb" sends left Boundry data to left  Neighbor "jb-1"
if (jb > 1) then     !99
   call MPI_send(u(its,jts-1), N, MPI_DOUBLE_PRECISION, Jb-1,&
&    LEFT, MPI_COMM_WORLD, IERR )
!STARTING AT ADDRESS u(ITS,JTS-1) FOR N ROWS.
!MPI step: taskid "jb" receives left Boundry data from left  Neighbor "jb-1"
if (jb > 1) then     !109
   call MPI_Recv(u(its,jts-1), N, &
 & MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RIGHT, MPI_COMM_WORLD, STATUS, ierr) 
!STARTING AT ADDRESS u(ITS,JTS-1) FOR N ROWS using 
!datatype mpi double precision, msgtag=RIGHT

!MPI step: taskid "jb" receives Right Boundry data from Right  Neighbor "jb+1"
if (jb<N ) then   !119
  call MPI_Recv(u(its,jts+1), N, &
 & MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, LEFT, MPI_COMM_WORLD, STATUS, ierr)

!STARTING AT ADDRESS u(ITS,JTS+1) FOR N ROWS using 
!datatype mpi double precision, msgtag=LEFT.

!MPI step: taskid "ib" sends TOP Boundry data to top Neighbor "ib-1"
if (ib<N) then     !129
   call MPI_send(u(its-1, jte), &
& N, MPI_DOUBLE_PRECISION, ib-1,TOP, MPI_COMM_WORLD, IERR )
!starting at address u(its-1,jte) for N columns.
!MPI step: taskid "ib" sends Bottom  Boundry data to Bottom  Neighbor "ib+1"
if (ib > 1) then     !139
   call MPI_send(u(its+1,jts), N, MPI_DOUBLE_PRECISION, ib+1,&
&    BOTTOM, MPI_COMM_WORLD, IERR )
!STARTING AT ADDRESS u(ITS+1,JTS) FOR N columns.

!MPI step: taskid "ib" receives left Boundry data from left  Neighbor "ib+1"
if (ib > 1) then     !149
   call MPI_Recv(u(its+1,jts), N, &
 & MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,TOP, MPI_COMM_WORLD, STATUS, ierr)
!STARTING AT ADDRESS u(ITS+1,JTS) FOR N columns using 
!datatype mpi double precision, msgtag=TOP
!MPI step: taskid "jb" receives Right Boundry data from Right  Neighbor "jb+1"
if (ib<N ) then   !159
  call MPI_Recv(u(its-1,jts), N, &
 & MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, BOTTOM, MPI_COMM_WORLD, STATUS, ierr)

!STARTING AT ADDRESS u(ITS-1,JTs) FOR N columnsS using 
!datatype mpi double precision, msgtag=BOTTOM.








!MPI STEP: task "jb" sends to MASTER LocalMaxAbsChange=err
! using mpi-max function with result stored at GlobalMaxAbsChange=gerr
!of count 1 and datatype mpi_double_precision.
  call MPI_Reduce(err, gerr,1,MPI_DOUBLE_PRECISION,&
& MPI_MAX,MASTER, MPI_COMM_WORLD, IERR )

!MASTER PRINT values every 100 iterations, 
!and check stopping Tolerance

if (taskid.eq. 0) then    !169
   if (mod(iter,100).eq. 0) then   !179
      write(10,1) taskid, iter, gerr
!      write(10,2) taskid, iter,ims,ime,jms,jme
 !print*, 'Exact Solution'
do i=ims,ime
    !write(10,3) 
   do j=jms,jme

end do
end do
     

 end if   !179
        if (gerr< tol) then  !189
           nstop=1
    end if      !189
end if      !169
!Broadcast the stopping value
  call MPI_Bcast(nstop,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,IERR)
1 format('taskid =', I4, ';iter= ', i5, ';GlobalMaxAbschange:gerr= ',e10.3) 
!3 format(('taskid=', I4, ';iter=', i5, ';f(' ,i3,',',i3,')=',3f9.4)


end if   !129
end if   !139 
end if   !149
end if   !159
end if    !119

end if      !109
 

end if     !99




end if    !88



 
End  Do !77


      call MPI_FINALIZE(ierr)










      end
!----------------------------------------------------------------
!-----------------------------------------------------------------
integer function divideup(m,n)  
! return m/n rounded up

integer :: ans 
divideup = ceiling(real(m)/real(n))
ans = (m+n-1)/n
! m =n : (m + n -1)/n = (2*n-1)/n = 1
! m =n+1 : (m + n -1)/n = (n+1 + n-1)/n = (2n)/n = 2
divideup = ans
return 
end function




