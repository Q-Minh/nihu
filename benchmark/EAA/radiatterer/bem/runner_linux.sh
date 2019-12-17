#! /bin/bash
export FREQ=$2
export FREQSTR=$FREQ"Hz"
export LE="100mm"
export TYPE=$1
export OUTDIR="data_fmm/"$TYPE"_"$LE
export EXE="./helmholtz_3d_hffmm_"$TYPE

$EXE --solve \
	--frequency $FREQ  \
	--surface_mesh "data/radiatterer_"$LE"_quad.off" \
	--surface_result $OUTDIR"/"$TYPE"_"$LE"_"$FREQSTR"_ps.res" \
	--cluster_size_constraint 0.2 
     
$EXE --postprocess --frequency $FREQ \
	 --surface_mesh "data/radiatterer_"$LE"_quad.off" \
	--field_mesh "data/radiatterer_points_quad.off" \
	--surface_result $OUTDIR"/"$TYPE"_"$LE"_"$FREQSTR"_ps.res" \
	--field_result $OUTDIR"/"$TYPE"_"$LE"_"$FREQSTR"_pf.res" \
	--cluster_size_constraint 0.2 
	
