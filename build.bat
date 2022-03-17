@echo off

:: For building the reference
echo +++ BUILDING REFERENCE...
cd reference && cargo build --release && echo +++ DONE BUILDING REFERENCE
cd ..
echo +++

:: For building the alignr
echo +++ BUILDING ALIGNR...
cd alignr && cargo build --release && echo +++ DONE BUILDING ALIGNR
cd ..
echo +++

:: For creating bin directory
echo +++ CREATING BIN DIRECTORY...
if exist bin (echo +++ `bin` directory already exists) else (md bin)

:: For placing the executbles in the bin directory
echo +++ COPYING REFERENCE EXECUTBLE TO BIN DIRECTORY
if exist bin\*.exe (del bin\*.exe && copy reference\target\release\ref-align.exe bin) else (copy reference\target\release\ref-align.exe bin)
echo +++ COPYING ALIGNR EXECUTBLE TO BIN DIRECTORY
copy alignr\target\release\alignr.exe bin
echo +++ DONE INITIALIZING BIN DIRECTORY
echo +++
