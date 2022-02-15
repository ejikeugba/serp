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

### Explanations on the Notes
- The suggested misspelled words in DESCRIPTION are all correct.
- The detritus check seems to be a false positive, no such directory or file exist in the temp directory. This appears to be common on rhub windows, especially with the recent Windows Server 2022.    



## win-builder

### Test environments
- using platform: x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2022-02-11 r81718 ucrt)

R CMD check results

Status: OK



## Further comments
- All dictated errors and warnings in the previous version are all addressed in the latest version of the package.

## Downstream dependencies
- There are currently no downstream dependencies for this package
