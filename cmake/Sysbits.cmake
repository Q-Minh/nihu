# Retrieve the number of bits (64 / 32) of opsystem
function(get_system_bits SYS_BITS)

# The 64-bit case
if(CMAKE_SIZEOF_VOID_P EQUAL 8) 
	set(${SYS_BITS} "64" PARENT_SCOPE)
# Otherwise 32 bit
else(CMAKE_SIZEOF_VOID_P EQUAL 8)
	set(${SYS_BITS} "32" PARENT_SCOPE)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

endfunction(get_system_bits)