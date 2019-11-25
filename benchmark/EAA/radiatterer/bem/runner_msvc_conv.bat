@ECHO OFF
SET FREQ=300
SET FREQSTR=%FREQ%Hz
SET LE=100mm
SET TYPE=const
SET OUTDIR=data_conv\%TYPE%_%LE%
SET EXE=helmholtz_3d_coll_%TYPE%.exe
ECHO ON

%EXE% --solve --frequency %FREQ%^
	--surface_mesh mesh\radiatterer_%LE%_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_ps.res
     
%EXE% --postprocess --frequency %FREQ%^
	--surface_mesh mesh\radiatterer_%LE%_quad.off^
	--field_mesh mesh\radiatterer_points_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_ps.res^
	--field_result %OUTDIR%\%TYPE%_%LE%_%FREQSTR%_pf.res
	