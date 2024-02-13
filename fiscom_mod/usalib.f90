     program libr
     use precision, only: pr => dp
     implicit none
     real(pr)      :: x,y
     integer       :: i
     integer(8)    :: j
     
     j = fftw_plan_dft_r2c_1d(128, x, y,i)
     
     end program 
     
