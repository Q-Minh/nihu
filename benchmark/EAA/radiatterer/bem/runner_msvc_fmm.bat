REM CALL: runner_msvc_fmm.bat frequency Le cluster_size_constraint
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
	--cluster_size_constraint %3^
	--restart 3000
     
%EXE% --postprocess --frequency %FREQUENCY%^
	--surface_mesh mesh\radiatterer_%LESTR%_quad.off^
	--field_mesh mesh\radiatterer_points_quad.off^
	--surface_result %OUTDIR%\%TYPE%_%LESTR%_%FREQSTR%_ps.res^
	--cluster_size_constraint %3^
	--field_result %OUTDIR%\%TYPE%_%LESTR%_%FREQSTR%_pf.res
	