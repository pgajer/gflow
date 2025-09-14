# valgrind

## Use Linux for leak checks (Docker/VM/GitHub Actions/R-hub). 
Use rhub 

## Xcode Instruments → Leaks (GUI, very good on macOS).

## leaks + malloc logging:

export MallocStackLogging=1   # enables stack logging
R -q -e 'your_repro_here' &   # run in background
leaks $(pgrep -n R)           # inspect the running R process

Note: AddressSanitizer works on macOS, but LeakSanitizer is disabled on macOS. So use Instruments/leaks locally, and keep Valgrind in CI on Linux.

