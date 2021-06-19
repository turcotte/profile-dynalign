/*-----------------------------------------------------------------------------
	to use ASSERT for debuging, only need to
	1. put '#define DEBUG' at the top of the program
	2. include this header file in the program
	3. comment out #define DEBUG to disable debuging
-----------------------------------------------------------------------------*/

#ifdef DEBUG
	void _Assert(char *, unsigned);    // prototype

	#define ASSERTFILE(str)		\
		static char strAssertFile[] = str ;
	#define ASSERT(f)	\
		if(f)						\
			{}						\
		else						\
			_Assert( __FILE__,  __LINE__ )


	void _Assert(char *strFile, unsigned uLine ) {
		fflush(NULL) ;
		fprintf(stderr, "\nAssertion failed: %s, line %u\n", strFile, uLine) ;
		fflush(stderr) ;
		abort() ;
	}

#else

	#define ASSERTFILE(str)
	#define ASSERT(f)

#endif

