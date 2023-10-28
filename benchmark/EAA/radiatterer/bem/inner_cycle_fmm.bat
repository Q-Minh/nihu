@REM USAGE: RUN_PARALLEL_FMM.BAT FROM STEP TO LE ClusterSize
@REM WHERE FROM AND TO DENOTE THE DOUBLE OF THE DESIRED VALUES

@ECHO OFF
SETLOCAL EnableDelayedExpansion
FOR /L %%F IN (%1, %2, %3) DO (
	SET /a "DOUBLE_FREQ=%%F"
	SET /a "FREQ_FLOOR=!DOUBLE_FREQ!/2"
	SET /a "FREQ_REM=!DOUBLE_FREQ! %% 2"
	IF /I "!FREQ_REM!" EQU "1" (
		SET FREQ=!FREQ_FLOOR!.5
	) ELSE (
		SET FREQ=!FREQ_FLOOR!
	)
	ECHO calling at frequency: !FREQ! Hz
	CALL runner_msvc_fmm !FREQ! %4 %5
)
