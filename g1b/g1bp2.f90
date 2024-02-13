program g1bp2
    use precision, only: pr=>dp
    use, intrinsic :: iso_c_binding

    implicit none
    include "/usr/include/fftw3.f03"
    type(C_ptr)                               :: plan
    character(len=60)                         :: fmt='(7(e28.21,5x))', tmp, loc="p2/"

    integer,parameter                         :: n=1300
    integer                                   :: i, j, k, n_ft, bin(101)

    real(pr), parameter                       :: pi=4._pr*atan(1._pr)
    real(pr)                                  :: factor, xmin, xmax, rmin, rmax, r_map, x_map, ly
    real(pr), allocatable                     :: x(:,:), r(:), in(:)

    complex(pr), allocatable                  :: out(:) 

    !######################################
    !           VAL INI
    !######################################

    allocate(x(4,n),r(5))
    n_ft = n-300
    allocate(in(n_ft),out(n_ft/2+1))
    
    x(:,1) = (/0.06_pr,0.3_pr,0.6_pr,0.9_pr/)
    r = (/1.5_pr,3.3_pr,3.5_pr,3.55_pr,4._pr/)

    do j=1,size(r)
        write(tmp,*) j
        plan = fftw_plan_dft_r2c_1d(n_ft,in(:),out(:),FFTW_MEASURE)

        !######################################
        !           CALCULO DE X
        !######################################

        open(33,file=trim(adjustl(loc))//"caos_"//trim(adjustl(tmp))//".d",status="replace")
        write(33,"(a15,e28.20)") "#r=", r(j)
        write(33,"(4(a15,15x))") "#x0=0.06", "x0=0.3", "x0=0.6", "x0=0.9"
            do i=1,n-1
                x(:,i+1)=(1._pr-x(:,i))*x(:,i)*r(j)
                if (i<=500) then
                    write(33,fmt) x(1,i), x(2,i), x(3,i), x(4,i)
                endif
            enddo
        close(33)

        !######################################
        !           BIN
        !######################################

        in(:)=x(3,300:)

        if (j==5) then
            open(35,file=trim(adjustl(loc))//"caos_bin.d",status="replace")
            write(35,"(a15,5x,a15)") "#r=4", "x0=0.6"
                xmax=in(1)
                xmin=in(1)

                do i=2,n_ft
                    if (in(i)>xmax) then
                        xmax=in(i)
                    else if (in(i)<xmin) then
                        xmin=in(i)
                    endif
                enddo

                bin(:)=0
                factor=(xmax-xmin)/real(size(bin),pr)
                write(35,"(a15,5x,e28.20)") "#", factor
                do i=1,size(bin)
                    do k=1,n_ft
                        if ((in(k)<=(xmin+real(i,pr)*factor)).and.(in(k)>(xmin+real(i-1,pr)*factor))) then
                            bin(i) = bin(i)+1
                        endif
                    enddo
                enddo
                write(35,'(i5)') bin(:)
            close(35)
        endif

        !######################################
        !           FFT
        !######################################

        open(34,file=trim(adjustl(loc))//"caos_ft_"//trim(adjustl(tmp))//".d",status="replace")
        write(34,"(4(a15,15x))") "#1/x", "real(F(1/x))", "img(F(1/x))"
            call fftw_execute_dft_r2c(plan, in(:), out(:))

            factor = n_ft-1._pr
            factor = 1._pr/factor
            do i = -n_ft/2,n_ft/2
            if ( i < 0) then
                write(34,fmt) factor*real(i,pr), &
                            real( conjg( out(-i+1) )/real(n_ft,pr) ), &
                            aimag(conjg( out(-i+1) )/real(n_ft,pr) )
            else if (i == 0) then
                write(34,fmt) factor*real(i,pr), &
                            real(out(1)/real(n_ft,pr)), &
                            aimag(out(1)/real(n_ft,pr))
            else
                write(34,fmt) factor*real(i,pr), &
                            real( out(i+1)/real(n_ft,pr) ), &
                            aimag( out(i+1)/real(n_ft,pr) )
            endif
            enddo
        close(34)
        call fftw_destroy_plan(plan)
    enddo

    !######################################
    !           MAPEO
    !######################################

    open(35,file=trim(adjustl(loc))//"caos_map.d",status="replace")
    write(35,"(4(a15,15x))") "#r", "x*"

    open(36,file=trim(adjustl(loc))//"caos_ly.d",status="replace")
    write(36,"(4(a15,15x))") "#r", "lambda"
        factor=0.001_pr
        j=0
        rmin=2._pr
        rmax=4._pr

        do while ((rmin+real(j,pr)*factor)<=rmax)
            r_map=rmin+real(j,pr)*factor
            x_map=0.01
            ly=0._pr

            do i=1,1500
                x_map=(1._pr-x_map)*x_map*r_map

                if ((i>=300).and.(i<=600).and.(r_map>=3.4_pr)) then
                    write(35,fmt) r_map, x_map
                else if (i>=300) then
                    ly=ly + log(abs(r_map*(1._pr-2._pr*x_map)))
                    if (i==1500) then
                        write(36,fmt) r_map, (ly/real(i,pr))
                    endif
                endif

            enddo
            j = j+1
        enddo
    close(35)
    close(36)

    !######################################
    !           MAPEO ZOOM
    !######################################

    open(37,file=trim(adjustl(loc))//"caos_mapzoom.d",status="replace")
    write(37,"(4(a15,15x))") "#r", "x*"
        factor=0.00001_pr
        j=0
        rmin=3.847_pr
        rmax=3.8568_pr

        do while ((rmin+real(j,pr)*factor)<=rmax)
            r_map=rmin+real(j,pr)*factor
            x_map=0.01
            xmin=0.44_pr
            xmax=0.56_pr

            do i=1,600
                x_map=(1-x_map)*x_map*r_map
                if ((i>300).and.(x_map<xmax).and.(x_map>xmin)) then
                    write(37,fmt) r_map, x_map
                endif
            enddo

            j = j+1
        enddo
    close(37)

end program g1bp2