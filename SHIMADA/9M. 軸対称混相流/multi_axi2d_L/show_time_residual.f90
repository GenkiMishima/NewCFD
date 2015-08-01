		subroutine show_time_residual
!
		include 'common.h'
		real rtime, ta(2)
!
			rtime=etime(ta)
			end_time=ta(1)
			step_cpu = end_time-start_time
			ellapse_time = end_time
			if(mod(ic,ipout).eq.0) then
			write(6,'(''ic='',i5,'' time='',e10.3,'' res='',e10.3,'' eCPU='',e10.3,'' cstep='',e10.3)') ic, &
						time, &
						log10(residual+1.e-20), &
						ellapse_time, &
						step_cpu
			end if

			write(70,'(i6,4e16.8)')  ic,time,log10(residual+1.e-20),ellapse_time,step_cpu
!
		return
		end