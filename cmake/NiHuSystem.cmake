# NiHuSystem.cmake 
# System related checks

# Variables set by this script
#   NIHU_SYS_BITS    System bits (32 or 64)

# Retrieve the number of bits (64 / 32) of operating system
function(get_system_bits SYS_BITS)
    # The 64-bit case
    if(CMAKE_SIZEOF_VOID_P EQUAL 8) 
        set(${SYS_BITS} "64" PARENT_SCOPE)
    # Otherwise 32 bit
    else()
        set(${SYS_BITS} "32" PARENT_SCOPE)
    endif()
endfunction(get_system_bits)
# end of function

message(STATUS "Checking system ...")
get_system_bits(NIHU_SYS_BITS)
message(STATUS "\t${NIHU_SYS_BITS}-bit system detected")
