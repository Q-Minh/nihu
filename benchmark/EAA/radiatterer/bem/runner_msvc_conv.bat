@ECHO OFF
SET FREQUENCY=%1
SET FREQSTR=%FREQUENCY%Hz
SET LE=100mm
SET TYPE=gauss
SET OUTDIR=data_conv\%TYPE%_%LE%
SET EXE=helmholtz_3d_coll_%TYPE%.exe
ECHO ON

%EXE% --solve --frequency %FREQUENCY%^
	--surface_mesh mesh\radiatterer_%LE%_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_ps.res^
	--restart 1000
     
%EXE% --postprocess --frequency %FREQUENCY%^
	--surface_mesh mesh\radiatterer_%LE%_quad.off^
	--field_mesh mesh\radiatterer_points_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_ps.res^
	--field_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_pf.res
