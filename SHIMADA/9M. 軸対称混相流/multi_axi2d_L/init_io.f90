		subroutine init_io
		include 'common.h'
!
		open(10,file='input.data',status='old')
		open(20,file='grid.data',status='old')
		open(50,file='field_bin.data',form='unformatted',status='unknown')
		open(51,file='density.data',status='unknown')
		open(52,file='temperature.data',status='unknown')
!		open(53,file='v.data',status='unknown')
		open(54,file='p_density.data',status='unknown')
		open(55,file='temp_p.data',status='unknown')
		open(60,file='input_echo.data',status='unknown')
!
		return
		end