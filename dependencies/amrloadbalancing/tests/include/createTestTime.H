#ifndef ROOT_CASEPATH
    #define ROOT_CASEPATH fileName("")
#endif

#ifndef ROOT_CASENAME
    #define ROOT_CASENAME fileName("")
#endif

#ifndef TIME_NAME
    #define TIME_NAME runTime
#endif

// Create Sample Time object with no function objects
Foam::Time TIME_NAME
(
        Foam::Time::controlDictName,
        ROOT_CASEPATH,
        ROOT_CASENAME,
        "system",
        "constant",
        false
);
