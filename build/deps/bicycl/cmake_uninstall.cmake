if(NOT EXISTS "/Users/tangguofeng/Documents/BBS-based DID/build/install_manifest.txt")
  message(FATAL_ERROR "Cannot find install manifest: /Users/tangguofeng/Documents/BBS-based DID/build/install_manifest.txt")
endif()

file(READ "/Users/tangguofeng/Documents/BBS-based DID/build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
foreach(file ${files})
  message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
  if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    execute_process(
      COMMAND "/opt/homebrew/Cellar/cmake/3.30.3/bin/cmake" -E rm "$ENV{DESTDIR}${file}"
      OUTPUT_VARIABLE rm_out
      RESULT_VARIABLE rm_retval
      )
    if(NOT "${rm_retval}" STREQUAL 0)
      message(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
    endif()
  else()
    message(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  endif()
endforeach()
