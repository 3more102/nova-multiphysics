# Validates target source lists at configure time. A duplicate can arise when
# independent feature branches each register the same implementation.
function(nova_assert_unique_target_sources target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR "NOVA source guard: unknown target '${target}'")
    endif()
    get_target_property(_sources "${target}" SOURCES)
    if(NOT _sources)
        return()
    endif()
    set(_seen "")
    foreach(_source IN LISTS _sources)
        if(_source IN_LIST _seen)
            message(FATAL_ERROR
                "NOVA source guard: '${target}' registers '${_source}' more than once")
        endif()
        list(APPEND _seen "${_source}")
    endforeach()
endfunction()
