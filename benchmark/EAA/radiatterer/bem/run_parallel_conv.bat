FOR /L %%S IN (1, 1, 12) DO (
	start inner_cycle_conv.bat %%S 12 2000
)
