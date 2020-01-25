    program po2d



    implicit none
    integer, allocatable :: psi(:,:), omg(:,:) ,tmp(:,:)
    integer                 i, j, k
    integer                 Nmesh, Stop_itr
    real                    Re, dt, eps


    write(*,*)'============== Grid generation ============='
    write(*,*) 'input meshsize'
    read(*,*) Nmesh

    allocate(psi(Nmesh,Nmesh))
    allocate(omg(Nmesh,Nmesh))
    allocate(tmp(Nmesh,Nmesh))

    psi(:,:)=0.0
    omg(:,:)=0.0
    tmp(:,:)=0.0

    write(*,*)'Mesh size=',Nmesh,'*',Nmesh,'=',Nmesh**2

    write(*,*)'============== Grid generation: Done ============='



    !MAIN LOOP
    do i = 1, Stop_itr

      !boundary condition for






      !close MAIN LOOP
    end do

    stop
    end
