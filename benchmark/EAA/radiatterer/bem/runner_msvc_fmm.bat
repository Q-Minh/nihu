@ECHO OFF
SET FREQUENCY=%1
SET FREQSTR=%FREQUENCY%Hz
SET LE=%2
SET LESTR=%LE%mm
SET TYPE=gauss
SET OUTDIR=data_fmm\%TYPE%_%LESTR%
SET EXE=helmholtz_3d_hffmm_%TYPE%.exe
ECHO ON

%EXE% --solve --frequency %FREQUENCY%^
	--surface_mesh mesh\radiatterer_%LESTR%_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LESTR%_%FREQSTR%_ps.res^
	--restart 1000
     
%EXE% --postprocess --frequency %FREQUENCY%^
	--surface_mesh mesh\radiatterer_%LESTR%_quad.off^
	--field_mesh mesh\radiatterer_points_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LESTR%_%FREQSTR%_ps.res^
	--field_result %OUTDIR%\%TYPE%_%LESTR%_%FREQSTR%_pf.res
	