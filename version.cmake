# check whether repository version can be determined
set(GIT_FOUND false)
if (UNIX)
	find_package(Git)
endif()

# get version variable via git or set it to "no-git-found"-value
if(GIT_FOUND)
	EXECUTE_PROCESS(
		COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --long --always
		OUTPUT_VARIABLE VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else()
	EXECUTE_PROCESS(
		COMMAND echo UNIX_AND_GIT_REQUIRED
		OUTPUT_VARIABLE VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
endif()



