REM helmholtz_3d_hf_fmm_standalone_msvc.exe^
	REM --surface_mesh data\radiatterer_050mm_quad.off^
	REM --surface_result data\const_quad_050mm\const_quad_050mm_500Hz_ps.res^
	REM --frequency 500^
	REM --solve

helmholtz_3d_hf_fmm_standalone_msvc.exe^
	--surface_mesh data\radiatterer_010mm_quad.off^
	--field_mesh data\radi_plane_010mm_quad.off^
	--surface_result data\const_quad_010mm\const_quad_010mm_4000Hz_ps.res^
	--field_result data\const_quad_010mm\const_quad_010mm_4000Hz_pf.res^
	--frequency 4000^
	--postprocess
