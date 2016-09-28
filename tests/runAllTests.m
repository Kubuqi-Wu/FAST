clc ; close all ; clear all ;

purgeTest = matlab.unittest.TestSuite.fromClass(?PurgeCutsTests) ;
run(purgeTest) ;

clear all ;

eVTest = matlab.unittest.TestSuite.fromClass(?EVTest) ;
run(eVTest) ;

clear all ;

mexTests = matlab.unittest.TestSuite.fromClass(?MexTests);
run(mexTests);

clear all ;

modelingTests = matlab.unittest.TestSuite.fromClass(?ModelingTests);
run(modelingTests);

clear all ;

paramsTests = matlab.unittest.TestSuite.fromClass(?ParamsTests) ;
run(paramsTests) ;

clear all ;

latticeTests = matlab.unittest.TestSuite.fromClass(?LatticeTests) ;
run(latticeTests) ;

clear all ;

compileLatticeTests = matlab.unittest.TestSuite.fromClass(?CompileLatticeTests);
run(compileLatticeTests);

clear all ;

sddpDeterministicTests = matlab.unittest.TestSuite.fromClass(?SddpDeterministicTest);
run(sddpDeterministicTests);

clear all ;

precutTests = matlab.unittest.TestSuite.fromClass(?PrecutTests);
run(precutTests);

clear all ;

dualsTests = matlab.unittest.TestSuite.fromClass(?DualsTests);
run(dualsTests);




