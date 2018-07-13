! INTERPOLATION OF NCEP FORCING FILES
! FORTRAN CODE BY RIESNA R AUDH

program interpolation

! This program interpolates the atmospheric forcing files
	implicit none

open(unit = 1, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\Cl.txt") !total cloud cover in percentage
open(unit = 2, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\Fsd_cloud.txt") !downward sw
open(unit = 3, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\P_rate.txt")       !(precipitation rate)
open(unit = 4, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\qa.txt")           !(specific humidity at 2 m height)
open(unit = 5, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\qs.txt")           !(specific humidity at the surface)
open(unit = 7, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\Ta.txt")           !(air temperature)
open(unit = 8, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\U.txt")            !(U wind speed at 10 m height)
open(unit = 9, file= "\\home\\riesna\\MSc2018\\ESIM_AntTest\\NCEP\\3\\V.txt")            !(V wind speed at 10 m height)


