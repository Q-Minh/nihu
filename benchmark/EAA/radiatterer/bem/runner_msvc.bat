@ECHO OFF
SET FREQ=600
SET FREQSTR=%FREQ%Hz
SET LE=050mm
SET TYPE=gauss
SET OUTDIR=data\const_quad_%LE%
SET EXE=helmholtz_3d_hffmm_%TYPE%.exe
ECHO ON

%EXE% --solve --frequency %FREQ%^
	--surface_mesh data\radiatterer_%LE%_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_ps.res
     
%EXE% --postprocess --frequency %FREQ%^
	--surface_mesh data\radiatterer_%LE%_quad.off^
	--field_mesh data\radi_plane_%LE%_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_ps.res^
	--field_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_pf.res
	