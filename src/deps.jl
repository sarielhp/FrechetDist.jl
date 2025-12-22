#! /bin/env julial

import Pkg;

Pkg.activate( "." )
Pkg.add( "TimerOutputs" )
Pkg.add( "Printf" )
Pkg.add( "CSV" )
Pkg.add( "DataFrames" )
Pkg.add( "PrettyTables" )
#Pkg.add( "CpuId" )
#Pkg.add( "Sys" )
Pkg.add( "CPUSummary" )
Pkg.add( "Hwloc" )
