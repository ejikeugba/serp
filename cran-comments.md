## Local CMD checks

### Test environments
- local x86_64-w64-mingw32 (64-bit), 4.1.2 (2021-11-01) 

R CMD check results
0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded


## R-hub builder

### Test environments
- R-hub Windows Server 2022, R-devel, 64 bit      (0 errors √ | 0 warnings √ | 2 notes x)
- R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC  (0 errors √ | 0 warnings √ | 1 notes x)
- R-hub Fedora Linux, R-devel, clang, gfortran    (0 errors √ | 0 warnings √ | 1 notes x)

R CMD check results

> On all platforms
    Possibly misspelled words in DESCRIPTION:
    CLM (3:47)
    SERP (12:53)
    Ugba (25:28, 26:5)
    al (26:13)
    et (26:10)
    
> On Windows Server 2022
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'


## win-builder

### Test environments
- using platform: x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2022-02-14 r81736 ucrt)

R CMD check results

> Possibly misspelled words in DESCRIPTION:
    CLM (3:47)
    SERP (12:53)
    Ugba (25:28, 26:5)
    al (26:13)
    et (26:10)

Status: 1 NOTE



## Downstream dependencies
- One reverse suggest having everything running smoothly


## Explanations on the Notes
- Suggested misspelled words in DESCRIPTION are all correct.
- The detritus check seems to be a false positive, no such directory or file exist in the temp directory. This appears to be common on rhub windows, especially with the recent Windows Server 2022.    



## Issues on previous version
- Reported errors and warnings in the previous version (v0.1.3), via continuous CRAN checks, have all been fixed in the latest version.

