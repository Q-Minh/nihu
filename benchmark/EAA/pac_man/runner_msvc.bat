REM CALL: runner_msvc.bat surface_mesh surface_excitation surface_result field_mesh field_result

@ECHO OFF
SET EXE=helmholtz_2d_wbfmm.exe
ECHO ON

%EXE% --solve^
	--surface_mesh %1^
	--surface_excitation %2^
	--surface_result %3^
	--num_leaf_nodes 20^
	--restart 3000
     
%EXE% --postprocess^
	--surface_mesh %1^
	--surface_excitation %2^
	--surface_result %3^
	--field_mesh %4^
	--field_result %5^
	--num_leaf_nodes 20
	